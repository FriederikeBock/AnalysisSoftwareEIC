#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

//*************************************************************************************
// Remove nonsense values at the end of the graph
//*************************************************************************************
void shrinkGraphErr(TGraphErrors* &graph, float energy){
  graph->Print();
  while(graph->GetN()>0 && graph->GetY()[graph->GetN()-1] == -1) graph->RemovePoint(graph->GetN()-1);
  while(graph->GetN()>0 && graph->GetY()[graph->GetN()-1] <= 0) graph->RemovePoint(graph->GetN()-1);
  while(graph->GetN()>0 && graph->GetX()[graph->GetN()-1] > energy) graph->RemovePoint(graph->GetN()-1);
  while(graph->GetN()>0 && graph->GetY()[0] == -1) graph->RemovePoint(0);
  while(graph->GetN()>0 && graph->GetY()[0] <= 0) graph->RemovePoint(0);
}

//*************************************************************************************
//*************************************************************************************
void ScaleByBinWidth2D(TH2F* h2){
  // Normalize by y bin size
  for(Int_t j = 1; j <= h2->GetNbinsX(); j++){
    for(Int_t i = 1; i <=  h2->GetNbinsY(); i++) {
      //printf("NLM2: i %d, j %d;  content %f / scale %f = %f\n",i,j,hNLM2->GetBinContent(j,i),scale2,hNLM2->GetBinContent(j,i)/scale2);
        h2->SetBinContent(j,i, h2->GetBinContent(j,i)/h2->GetXaxis()->GetBinWidth(j));
        h2->SetBinError  (j,i, h2->GetBinError  (j,i)/h2->GetXaxis()->GetBinWidth(j));
    } // y bin loop 
  } // x bin loop
  
}


//*************************************************************************************
//*************************************************************************************
void pidreso_Pythia(
                            TString inputFileName   = "file.root",
                            TString suffix          = "pdf",
                            TString addLabel        = "",
                            Double_t tRes           = 25.,          
                            Bool_t properFit        = kTRUE,
                            Int_t trackCuts         = 0
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  SetPlotStyle();

  TString dateForOutput             = ReturnDateStringForOutput();
  TString readTrackClass            = "";
  TString readRegion[3]             = {"ETTL", "CTTL", "FTTL"};
  TString readTimeDet[2]            = {"wScatE", "woScatE"};
  TString readPIDSpecial[3]         = {"pi", "K", "p"};
  TString labelPIDSpecial[3]         = {"#pi^{#pm}", "K^{#pm}", "p/#bar{p}"};
  
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
  TString outputDirOverview         = Form("plots/%s/%s/Overview",dateForOutput.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDirBetaRes);
  gSystem->Exec("mkdir -p "+outputDirBeta);
  gSystem->Exec("mkdir -p "+outputDirOverview);
  
  TString perfLabel = "#it{#bf{ECCE}} simulation";
  TString collisionSystem = "Pythia 6, e+p, 18x275 GeV";
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addLabel.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  } else if (addLabel.Contains("Q2100")) {
    pTHard = "#it{Q}^{2} > 100 GeV^{2}";
    nLinesCol++;
  }

  Int_t nActiveEta            = 14;
  Double_t maxBetaSigma       = 0.05;
  if (addLabel.Contains("Hist")){
    maxBetaSigma         = 0.05;
  }
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;


  //************************** Read data **************************************************
  TH1D* h_mean_p_betaReso[6][nEta+1]            = {{NULL}};
  TH1D* h_sigma_p_betaReso[6][nEta+1]           = {{NULL}};  
  TH1D* h_beta_reso_bins[6][nEta+1][nPLow]        = {{{NULL}}};
  TF1* fit_tracks_reso_bins[6][nEta+1][nPLow]     = {{{NULL}}};

  TH1D* h_mean_p_betaSmearT0Reso[6][nEta+1]        = {{NULL}};
  TH1D* h_sigma_p_betaSmearT0Reso[6][nEta+1]       = {{NULL}};  
  TH1D* h_betaSmearT0_reso_bins[6][nEta+1][nPLow]    = {{{NULL}}};
  TF1* fit_tracksSmearT0_reso_bins[6][nEta+1][nPLow] = {{{NULL}}};
  
  TH2F* h_beta_reso[6][nEta+1]            = {{NULL}};
  TH2F* h_betaSmearT0_reso[6][nEta+1]        = {{NULL}};
  TH2F* h_beta_p[6][nEta+1]               = {{NULL}};
  TH2F* h_betaGen_p[6][nEta+1]            = {{NULL}};
  TH2F* h_betaTrack_p[6][nEta+1]          = {{NULL}};
  TH2F* h_betaSmearT0_p[6][nEta+1]          = {{NULL}};
  
  TH2F* h_beta_p_Region[6][3]             = {{NULL}};
  TH2F* h_betaSmearT0_p_Region[6][3]      = {{NULL}};
  TH2F* h_NSigma_p_Region[6][6][3]        = {{NULL}};   // [currentpart][refpart][region] 
  TGraphErrors* g_ParticleYieldTot_p_Region[6][6][3]      = {NULL}; // [currentpart][refpart][region] 
  TGraphErrors* g_ParticleYieldCut_p_Region[6][6][3]      = {NULL}; // [currentpart][refpart][region] 
  TGraphErrors* g_ParticleRejection_p_Region[6][6][3]      = {NULL}; // [currentpart][refpart][region] 
  
  
  TH2F* h_betaTrack_p_Region[6][3]        = {{NULL}};
  TH2F* h_betaGen_p_Region[6][3]          = {{NULL}};
  
  TH2F* h_InitT0VsNtrk[2]                 = {NULL};
  TH2F* h_T0VsNtrk[2]                     = {NULL};
  TH2F* h_InitT0VsNttl[2]                 = {NULL};
  TH2F* h_T0VsNttl[2]                     = {NULL};
  TH2F* h_InitT0DiffVsNtrk[2]             = {NULL};
  TH2F* h_T0DiffVsNtrk[2]                 = {NULL};
  TH2F* h_InitT0DiffVsNttl[2]             = {NULL};
  TH2F* h_T0DiffVsNttl[2]                 = {NULL};
  TH2F* h_NttlvsNtrk_initT0[2]            = {NULL};
  TH2F* h_NttlvsNtrk_T0[2]                = {NULL};
  
  TH2F* h_InitDbetaVsP[2][3]              = {{NULL}};
  TH2F* h_DbetaVsP[2][3]                  = {{NULL}};
  
  TH2F* hBetaVsPDet[3]                    = {NULL};
  TH2F* hResBetaVsPDet_wT0[3][6]          = {NULL};
  TH2F* hResBetaVsPDet_woT0[3][6]         = {NULL};
  TGraphErrors* g_betaMeanSmearT0_p_Region[6][3]      = {NULL};
  TGraphErrors* g_betaSigmaAbsSmearT0_p_Region[6][3]  = {NULL};
  TGraphErrors* g_betaSigmaDownSmearT0_p_Region[6][3] = {NULL};
  TGraphErrors* g_betaSigmaUpSmearT0_p_Region[6][3]   = {NULL};
  
  
  Int_t nSlicesP                          = 5;
  Float_t sliceP[5]                       = {0.3, 1.5, 2., 3., 4.5};
  TH1D* hSliceBetaRegion[6][3][5]         = {{{NULL}}};
  
  bool enablePIDEta[6][nEta+1]            = {{0}};
  bool enablePIDEtaSmear[6][nEta+1]       = {{0}};
  Int_t nEvents                           = 0;
  Color_t colorIter[4]                    = {kGray+2, kBlue+1, kGreen+2, kRed+1};
  Style_t styleMarkerIter[4]              = {24, 27, 28, 20};
  Size_t sizeMarkerIter[4]                = {1.2, 1.7, 1.2, 1.3};
    
  TFile* inputFile  = new TFile(inputFileName.Data());
  
  TH2D* h_T0VsIter              = (TH2D*)inputFile->Get("h_T0DiffVsIter");
  TH1D* h_NEvents               = (TH1D*)inputFile->Get("h_TPH_NEvents");
  nEvents = h_NEvents->GetBinContent(1);
  TH2D* h_UPartIter             = (TH2D*)inputFile->Get("h_TPH_UsedParticlesVsIter");
  Int_t nIter                   = 4;
  TH1D* h_T0VsIterProj[2*nIter];
  TH1D* h_UPartIterProj[2*nIter];
  for (Int_t it = 0; it < 2*nIter; it++){
    h_T0VsIterProj[it]          = (TH1D*)h_T0VsIter->ProjectionY(Form("h_V0VsIter_%i", it), it+1, it+1 ,"E");
    h_UPartIterProj[it]         = (TH1D*)h_UPartIter->ProjectionY(Form("h_UPartIter_%i", it), it+1, it+1 ,"E");
  }

  for (Int_t td = 0; td < 2; td++){
    h_InitT0VsNtrk[td]                    = (TH2F*)inputFile->Get(Form("h_InitT0VsNtrk_%s", readTimeDet[td].Data()));
    h_T0VsNtrk[td]                        = (TH2F*)inputFile->Get(Form("h_T0VsNtrk_%s", readTimeDet[td].Data()));
    h_InitT0VsNttl[td]                    = (TH2F*)inputFile->Get(Form("h_InitT0VsNttlhit_%s", readTimeDet[td].Data()));
    h_T0VsNttl[td]                        = (TH2F*)inputFile->Get(Form("h_T0VsNttlhit_%s", readTimeDet[td].Data()));
    h_InitT0DiffVsNtrk[td]                = (TH2F*)inputFile->Get(Form("h_InitT0DiffVsNtrk_%s", readTimeDet[td].Data()));
    h_T0DiffVsNtrk[td]                    = (TH2F*)inputFile->Get(Form("h_T0DiffVsNtrk_%s", readTimeDet[td].Data()));
    h_InitT0DiffVsNttl[td]                = (TH2F*)inputFile->Get(Form("h_InitT0DiffVsNttlhit_%s", readTimeDet[td].Data()));
    h_T0DiffVsNttl[td]                    = (TH2F*)inputFile->Get(Form("h_T0DiffVsNttlhit_%s", readTimeDet[td].Data()));
    h_NttlvsNtrk_initT0[td]               = (TH2F*)inputFile->Get(Form("h_NttlhitVsNtrk_wInitT0_%s", readTimeDet[td].Data()));
    h_NttlvsNtrk_T0[td]                   = (TH2F*)inputFile->Get(Form("h_NttlhitVsNtrk_wFinalT0_%s", readTimeDet[td].Data()));
    for (Int_t pidS = 0; pidS < 3; pidS++){
      h_InitDbetaVsP[td][pidS]            = (TH2F*)inputFile->Get(Form("h_InitDbetaVsP_%s_%s", readPIDSpecial[pidS].Data(), readTimeDet[td].Data()));
      h_DbetaVsP[td][pidS]                = (TH2F*)inputFile->Get(Form("h_DbetaVsP_%s_%s", readPIDSpecial[pidS].Data(), readTimeDet[td].Data()));
    }
  }
  
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    for (Int_t pid = 0; pid < 6; pid++){    
      h_beta_reso[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_Res_InvSmearBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      if (!h_beta_reso[pid][iEta]){
        cout << "tried to find: " << Form("h_Res_InvSmearBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta) << " and failed" << endl;
        continue;
      } else {
        enablePIDEta[pid][iEta] = true;
//         cout <<  "found: " << Form("h_Res_InvSmearBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta) << endl;
      }
      h_beta_reso[pid][iEta]->Sumw2();
      h_mean_p_betaReso[pid][iEta] = new TH1D(Form("histBetaResol%s_%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                            ";#it{p}^{MC} (GeV/#it{c}); #LT (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} #GT",
                                            nPLow, partP);
      h_sigma_p_betaReso[pid][iEta] = new TH1D(Form("histBetaResol%s_%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                              ";#it{p}^{MC} (GeV/#it{c}); #sigma( (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} )", 
                                              nPLow, partP);
      h_betaSmearT0_reso[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_Res_InvSmearAndT0Beta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      if (h_betaSmearT0_reso[pid][iEta]){
        enablePIDEtaSmear[pid][iEta]   = true;
      }
      if (enablePIDEtaSmear[pid][iEta]){
        h_betaSmearT0_reso[pid][iEta]->Sumw2();
        h_mean_p_betaSmearT0Reso[pid][iEta] = new TH1D(Form("histBetaResol%sWT0%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                              ";#it{p}^{MC} (GeV/#it{c}); #LT (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} #GT",
                                              nPLow, partP);
        h_sigma_p_betaSmearT0Reso[pid][iEta] = new TH1D(Form("histBetaResol%sWT0%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta), 
                                                ";#it{p}^{MC} (GeV/#it{c}); #sigma( (1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC} )", 
                                                nPLow, partP);
        h_betaSmearT0_p[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvSmearAndT0Beta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
        h_betaSmearT0_p[pid][iEta]->Sumw2();
      }
      h_betaTrack_p[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_betaTrack_p[pid][iEta]->Sumw2();
      h_betaGen_p[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvGenBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_betaGen_p[pid][iEta]->Sumw2();
      h_beta_p[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_InvSmearBeta%s_Eta_p_%s_%d", readTrackClass.Data(), partName[pid].Data(), iEta));
      h_beta_p[pid][iEta]->Sumw2();
    }
  }
  
  for (Int_t eR = 0; eR < 3; eR++){
    hBetaVsPDet[eR]                     = (TH2F*)inputFile->Get(Form("h_BetaVsP_T0_%s", readRegion[eR].Data()));
    
    for (Int_t pid = 0; pid < 6; pid++){    
      hResBetaVsPDet_wT0[eR][pid]       = (TH2F*)inputFile->Get(Form("h_DbetaVsP_wT0_%s_%s", readRegion[eR].Data(), partName[pid].Data()));
      hResBetaVsPDet_woT0[eR][pid]      = (TH2F*)inputFile->Get(Form("h_DbetaVsP_woT0_%s_%s", readRegion[eR].Data(), partName[pid].Data()));
      
      cout << "region " << eR << "adding " ;
      for(Int_t iEta=minEtaBinTTL[eR]; iEta<maxEtaBinTTL[eR]+1;iEta++){
        if (!enablePIDEta[pid][iEta]) continue;
        cout << iEta << "\t";
        if (!h_beta_p_Region[pid][eR]){
          h_beta_p_Region[pid][eR]      = (TH2F*)h_beta_p[pid][iEta]->Clone(Form("h_InvSmearBeta%s_Eta_p_%s_%s", readTrackClass.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data()));
          h_beta_p_Region[pid][eR]->Sumw2();
          h_betaTrack_p_Region[pid][eR] = (TH2F*)h_betaTrack_p[pid][iEta]->Clone(Form("h_InvBeta%s_Eta_p_%s_%s", readTrackClass.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data()));
          h_betaTrack_p_Region[pid][eR]->Sumw2();
          h_betaGen_p_Region[pid][eR]   = (TH2F*)h_betaGen_p[pid][iEta]->Clone(Form("h_InvGenBeta%s_Eta_p_%s_%s", readTrackClass.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data()));
          h_betaGen_p_Region[pid][eR]->Sumw2();
        } else {
          h_beta_p_Region[pid][eR]->Add(h_beta_p[pid][iEta]);
          h_betaTrack_p_Region[pid][eR]->Add(h_betaTrack_p[pid][iEta]);
          h_betaGen_p_Region[pid][eR]->Add(h_betaGen_p[pid][iEta]);
        }
        
        if (enablePIDEtaSmear[pid][iEta]){
          if (!h_betaSmearT0_p_Region[pid][eR]){
            h_betaSmearT0_p_Region[pid][eR]      = (TH2F*)h_betaSmearT0_p[pid][iEta]->Clone(Form("h_InvSmearT0Beta%s_Eta_p_%s_%s", readTrackClass.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data()));
            h_betaSmearT0_p_Region[pid][eR]->Sumw2();
          } else {
            h_betaSmearT0_p_Region[pid][eR]->Add(h_betaSmearT0_p[pid][iEta]);
          }
        } 
      }
      cout << endl;
      cout << "eta range: " << partEta[minEtaBinTTL[eR]] << "\t" << partEta[maxEtaBinTTL[eR]+1] << endl;
      for (Int_t pt = 0; pt < nSlicesP; pt++){
        hSliceBetaRegion[pid][eR][pt] = (TH1D*)h_betaSmearT0_p_Region[pid][eR]->ProjectionY(Form("sliceP_%i_%i_%i", eR, pid, pt), h_betaSmearT0_p_Region[pid][eR]->GetXaxis()->FindBin(sliceP[pt]), h_betaSmearT0_p_Region[pid][eR]->GetXaxis()->FindBin(sliceP[pt]), "E");
      }
//       if (pid == 1){
//         TH2F* temp = (TH2F*)h_betaSmearT0_p_Region[pid][eR]->RebinX(4,Form("h_InvSmearT0Beta%s_Eta_p_%s_%s", readTrackClass.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data()));
//         h_betaSmearT0_p_Region[pid][eR] = temp;
//       }
      if (pid > 0){
        Double_t mean[500];
        Double_t sigma[500];
        Double_t sigmaU[500];
        Double_t sigmaD[500];
        Double_t p[500];
        Double_t pE[500];
        Double_t meanE[500];
        Double_t sigmaE[500];
        Int_t bins = 0;

        TH2F* temp = (TH2F*)h_betaSmearT0_p_Region[pid][eR]->Clone("toslice");
        if (pid == 1){
          temp = (TH2F*)h_betaSmearT0_p_Region[pid][eR]->RebinX(4,"toslice");
        }

        
        for (Int_t iPt = 1; iPt < temp->GetNbinsX() && bins < 500; iPt++ ){
          TH1D* tempProjec = (TH1D*)temp->ProjectionY("dummy", iPt, iPt, "E");
          if (pid == 4 && (p[bins] < 0.55 || p[bins] > 2))
            tempProjec->Rebin(4);
          if (pid == 5 &&  p[bins] > 2)
            tempProjec->Rebin(4);
          if (pid == 1 &&  p[bins] < 1)
            tempProjec->Rebin(4);
          if (pid == 3 &&  p[bins] < 0.2)
            tempProjec->Rebin(4);
          if (tempProjec->GetEntries() > 50){
            TF1* tempGauss = new TF1("tempGauss",  "gaus", 0.0905,1.495);
            tempGauss->SetParameter(1,tempProjec->GetMean());
            tempGauss->SetParameter(2,tempProjec->GetStdDev());
            tempProjec->Fit(tempGauss,"QNRME+","",tempProjec->GetMean()-2*tempProjec->GetStdDev(),tempProjec->GetMean()+2*tempProjec->GetStdDev());
            
            mean[bins]  = tempGauss->GetParameter(1);
            meanE[bins] = tempGauss->GetParError(1);
            p[bins]     = temp->GetXaxis()->GetBinCenter(iPt);
            pE[bins]    = temp->GetXaxis()->GetBinWidth(iPt)/2;
            sigma[bins]   = tempGauss->GetParameter(2);
            sigmaU[bins]  = mean[bins]+3*tempGauss->GetParameter(2);
            sigmaD[bins]  = mean[bins]-3*tempGauss->GetParameter(2);
            sigmaE[bins] = tempGauss->GetParError(2);
            if (pid == 3 && p[bins] < 0.14 ) continue;
            if (pid == 4 && p[bins] < 0.48 ) continue;
            if (pid == 5 && p[bins] < 0.85 ) continue;
            bins++;
          }
        }
        g_betaMeanSmearT0_p_Region[pid][eR]   = new TGraphErrors(bins, p, mean, pE, meanE);
        g_betaSigmaAbsSmearT0_p_Region[pid][eR]  = new TGraphErrors(bins, p, sigma, pE, sigmaE);
        g_betaSigmaUpSmearT0_p_Region[pid][eR]  = new TGraphErrors(bins, p, sigmaU, pE, sigmaE);
        g_betaSigmaDownSmearT0_p_Region[pid][eR]  = new TGraphErrors(bins, p, sigmaD, pE, sigmaE);
      }
    }
    
    cout << "here" << endl;
    for (Int_t pid = 0; pid < 6; pid++){
       for (Int_t pidref = 1; pidref < 6; pidref++){
         cout << __LINE__ <<"\t" << pid <<"\t"<< pidref<< endl;
          h_NSigma_p_Region[pid][pidref][eR] = new TH2F(Form("h_%s_NSigma%s_Eta_p_%s", partName[pid].Data(), partName[pidref].Data(), nameOutEtaRange[eR].Data()),"", nBinsPFine, binningPFine ,200, -15, 15);
          if (h_betaSmearT0_p_Region[pid][eR] && g_betaMeanSmearT0_p_Region[pidref][eR]){
            for (Int_t ip = 1; ip < h_betaSmearT0_p_Region[pid][eR]->GetNbinsX()+1; ip++){
              for (Int_t ibe = 0; ibe < h_betaSmearT0_p_Region[pid][eR]->GetNbinsY()+2; ibe++){
                Double_t p      = h_betaSmearT0_p_Region[pid][eR]->GetXaxis()->GetBinCenter(ip);
                Double_t nsigma = (h_betaSmearT0_p_Region[pid][eR]->GetYaxis()->GetBinCenter(ibe)-g_betaMeanSmearT0_p_Region[pidref][eR]->Eval(p))/g_betaSigmaAbsSmearT0_p_Region[pidref][eR]->Eval(p);
                h_NSigma_p_Region[pid][pidref][eR]->Fill(p, nsigma, h_betaSmearT0_p_Region[pid][eR]->GetBinContent(ip,ibe));
              }
            }
          }
          
          if (pid > 0 && pid != pidref && h_NSigma_p_Region[pid][pidref][eR]){
            Double_t yieldTot[500];
            Double_t yieldErrTot[500];
            Double_t yieldCut[500];
            Double_t yieldErrCut[500];
            Double_t reject[500];
            Double_t rejectErr[500];
            Double_t p[500];
            Double_t pE[500];
            Int_t bins = 0;
            
            for (Int_t iPt = 0; iPt < nBinsPLow ; iPt++ ){
              TH1D* tempProjec = (TH1D*)h_NSigma_p_Region[pid][pidref][eR]->ProjectionY("dummy", h_NSigma_p_Region[pid][pidref][eR]->GetXaxis()->FindBin(binningP[iPt]+0.0001), h_NSigma_p_Region[pid][pidref][eR]->GetXaxis()->FindBin(binningP[iPt+1]-0.0001), "E");
              
              
          
            yieldTot[bins]      = tempProjec->Integral();
            yieldErrTot[bins]   = TMath::Sqrt(yieldTot[bins]);
            yieldCut[bins]      = tempProjec->Integral(tempProjec->GetXaxis()->FindBin( -3+0.001), tempProjec->GetXaxis()->FindBin( 3+0.001));
            yieldErrCut[bins]   = TMath::Sqrt(yieldCut[bins]);
            if (yieldCut[bins] > 0){
              reject[bins]      = yieldTot[bins]/yieldCut[bins];
              rejectErr[bins]   = reject[bins]*TMath::Sqrt( TMath::Power(yieldErrTot[bins]/yieldTot[bins],2) + TMath::Power(yieldErrCut[bins]/yieldCut[bins],2));
            } else {
              reject[bins]      = -1;
              rejectErr[bins]   = 1;
            }
            p[bins]     = (binningP[iPt]+binningP[iPt+1])/2;
            pE[bins]    = (binningP[iPt+1]-binningP[iPt])/2;
            bins++;
          }
          g_ParticleYieldTot_p_Region[pid][pidref][eR]      = new TGraphErrors(bins, p, yieldTot, pE, yieldErrTot);
          g_ParticleYieldCut_p_Region[pid][pidref][eR]      = new TGraphErrors(bins, p, yieldCut, pE, yieldErrCut);
          g_ParticleRejection_p_Region[pid][pidref][eR]     = new TGraphErrors(bins, p, reject, pE, rejectErr);
          cout << "rejection " << pid << "\t" << partName[pid].Data() << "\t ref \t" << pidref  << "\t" << partName[pidref].Data() << endl;
//           cout << "tot" << endl;
//           g_ParticleYieldTot_p_Region[pid][pidref][eR]->Print();
//           cout << "cut" << endl;
//           g_ParticleYieldCut_p_Region[pid][pidref][eR]->Print();
//           cout << "rejection" << endl;
//           g_ParticleRejection_p_Region[pid][pidref][eR]->Print();
          
        }
      }
    }
  }

  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput)std::cout << Form("%1.1f < #eta < %1.1f",etaMin,etaMax)  << std::endl;
    for (Int_t iPt = 0; iPt < nPLow; iPt++){
      for (Int_t pid = 0; pid < 6; pid++){
        if (!enablePIDEta[pid][iEta]) continue;
  //       std::cout << partP[iPt] << "\t" << h_beta_reso[iEta]->GetXaxis()->FindBin(partP[iPt]) << std::endl;
        h_beta_reso_bins[pid][iEta][iPt] = (TH1D*)h_beta_reso[pid][iEta]->ProjectionY(Form("projection_%s_dBeta_%d_%d",partName[pid].Data(),iPt, iEta), 
                                                                          h_beta_reso[pid][iEta]->GetXaxis()->FindBin(partP[iPt]), h_beta_reso[pid][iEta]->GetXaxis()->FindBin(partP[iPt+1]),"e"); 
        h_beta_reso_bins[pid][iEta][iPt]->Rebin(rebinEta_BetaResol[iEta]);
        h_beta_reso_bins[pid][iEta][iPt]->GetXaxis()->SetRangeUser(-0.1,0.1);
        if(h_beta_reso_bins[pid][iEta][iPt]->GetEntries() > 10 && properFit){
            fit_tracks_reso_bins[pid][iEta][iPt] = new TF1(Form("fitGauss%s%d_%d",partName[pid].Data() ,iPt, iEta),  "gaus", -betaResolFit[iEta],betaResolFit[iEta]);
            fit_tracks_reso_bins[pid][iEta][iPt]->SetParameter(0,h_beta_reso_bins[pid][iEta][iPt]->GetMaximum() );
            fit_tracks_reso_bins[pid][iEta][iPt]->SetParLimits(0,0.9*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum(),4*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum() );
            h_beta_reso_bins[pid][iEta][iPt]->Fit(fit_tracks_reso_bins[pid][iEta][iPt],"QNRME+","",-betaResolFit[iEta],betaResolFit[iEta]);
        }
        
        if (enablePIDEtaSmear[pid][iEta]){  
          h_betaSmearT0_reso_bins[pid][iEta][iPt] = (TH1D*)h_betaSmearT0_reso[pid][iEta]->ProjectionY(Form("projectionAEMC_%s_dBeta_%d_%d",partName[pid].Data(),iPt, iEta), 
                                                                            h_betaSmearT0_reso[pid][iEta]->GetXaxis()->FindBin(partP[iPt]), h_betaSmearT0_reso[pid][iEta]->GetXaxis()->FindBin(partP[iPt+1]),"e"); 
          h_betaSmearT0_reso_bins[pid][iEta][iPt]->Rebin(rebinEta_BetaResol[iEta]);
          h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetXaxis()->SetRangeUser(-0.1,0.1);
          if(h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetEntries() > 10 && properFit){
              fit_tracksSmearT0_reso_bins[pid][iEta][iPt] = new TF1(Form("fitGaussAEMC%s%d_%d",partName[pid].Data() ,iPt, iEta),  "gaus", -betaResolFit[iEta],betaResolFit[iEta]);
              fit_tracksSmearT0_reso_bins[pid][iEta][iPt]->SetParameter(0,h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetMaximum() );
              fit_tracksSmearT0_reso_bins[pid][iEta][iPt]->SetParLimits(0,0.9*h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetMaximum(),4*h_beta_reso_bins[pid][iEta][iPt]->GetMaximum() );
              h_betaSmearT0_reso_bins[pid][iEta][iPt]->Fit(fit_tracksSmearT0_reso_bins[pid][iEta][iPt],"QNRME+","",-betaResolFit[iEta],betaResolFit[iEta]);
          }
        }
      }
      
      if(iEta == 0 && debugOutput){
        std::cout <<h_mean_p_betaReso[0][iEta]->GetBinCenter(iPt+1) << "\t"<<iPt+1<<  "\t"<< h_beta_reso_bins[0][iEta][iPt]->GetMean() <<  "\t"<< h_beta_reso_bins[0][iEta][iPt]->GetRMS() << std::endl;
        if (fit_tracks_reso_bins[0][iEta][iPt]) std::cout<< fit_tracks_reso_bins[0][iEta][iPt]->GetParameter(1) << "\t" << fit_tracks_reso_bins[0][iEta][iPt]->GetParameter(2) << std::endl;
      }
      if (properFit){
        for (Int_t pid = 0; pid < 6; pid++){
          if (!enablePIDEta[pid][iEta]) continue;
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
          if (enablePIDEtaSmear[pid][iEta]){  
            if (fit_tracksSmearT0_reso_bins[pid][iEta][iPt]){
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,fit_tracksSmearT0_reso_bins[pid][iEta][iPt]->GetParameter(1));
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,fit_tracksSmearT0_reso_bins[pid][iEta][iPt]->GetParError(1));
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,fit_tracksSmearT0_reso_bins[pid][iEta][iPt]->GetParameter(2));
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,fit_tracksSmearT0_reso_bins[pid][iEta][iPt]->GetParError(2));
            } else {
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,-100000);
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,0);
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,-1000);
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,0);         
            }
          }
        }
      } else {
        for (Int_t pid = 0; pid < 6; pid++){
          if (!enablePIDEta[pid][iEta]) continue;
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
          if (enablePIDEtaSmear[pid][iEta]){  
            if(h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetEntries() > 10){
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetMean());
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetMeanError());
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetRMS());
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetRMSError());      
  //       h_sigma_p_betaReso[iEta]->SetBinContent(iPt+1,h_beta_reso_bins[iEta][iPt]->GetStdDev());
  //       h_sigma_p_betaReso[iEta]->SetBinError(iPt+1,h_beta_reso_bins[iEta][iPt]->GetStdDevError());
            } else {
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,-100000);
              h_mean_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,0);
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinContent(iPt+1,-1000);
              h_sigma_p_betaSmearT0Reso[pid][iEta]->SetBinError(iPt+1,0);         
            }
          }
        }  
      }
    }
  }

  cout << __LINE__ << endl;
  
  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  split_canvas(cPNG, "cPNG1", nPLow);
  
  TCanvas* cSingle2D = new TCanvas("cSingle2D","",0,0,1000,800);
  DrawGammaCanvasSettings( cSingle2D, 0.1, 0.12, 0.025, 0.105);
  cSingle2D->SetLogz();
  cSingle2D->SetLogx();
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Int_t padnum =0;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput) std::cout << Form("%1.1f < #eta < %1.1f",etaMin,etaMax)  << std::endl;
      
    for (Int_t pid = 0; pid < 6; pid++){
      if (!enablePIDEta[pid][iEta]) continue;
      padnum =0;
      for(Int_t iPt=0; iPt< nPLow; iPt++){
        cPNG->cd(padnum+1);
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1)->SetLogy(1);
//        cPNG->cd(padnum+1)->SetLogx(1);
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
          drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        
        drawLatexAdd(Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f ",partP[iPt],partP[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f \t",partP[iPt],partP[iPt+1]) << h_mean_p_betaReso[pid][iEta]->GetBinCenter(iPt+1)<< "\t"<< h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1) << std::endl;
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
      
      if (enablePIDEtaSmear[pid][iEta]){  
        padnum =0;
        for(Int_t iPt=0; iPt< nPLow; iPt++){
          cPNG->cd(padnum+1);
          cPNG->cd(padnum+1)->Clear();
          padnum++;
        }
        padnum =0;
        for(Int_t iPt=0; iPt< nPLow; iPt++){
          cPNG->cd(padnum+1);
          DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
          cPNG->cd(padnum+1)->SetLogy(1);
          SetStyleHistoTH1ForGraphs(h_betaSmearT0_reso_bins[pid][iEta][iPt], "(1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC}", "counts",
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
          h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetYaxis()->SetRangeUser(0.8,h_beta_reso_bins[pid][iEta][iPt]->GetMaximum()*2);
          
          DrawGammaSetMarker(h_betaSmearT0_reso_bins[pid][iEta][iPt], 20, 1.5, kBlack, kBlack);      
          h_betaSmearT0_reso_bins[pid][iEta][iPt]->Draw("p,e");
          
          if (fit_tracksSmearT0_reso_bins[pid][iEta][iPt] && properFit){
            DrawGammaSetMarkerTF1( fit_tracksSmearT0_reso_bins[pid][iEta][iPt], 7, 2, kBlue+1);
            fit_tracksSmearT0_reso_bins[pid][iEta][iPt]->Draw("same");
          }
          
          if (padnum == 0){
            drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
            drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          }
          
          drawLatexAdd(Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f ",partP[iPt],partP[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f \t",partP[iPt],partP[iPt+1]) << h_mean_p_betaReso[pid][iEta]->GetBinCenter(iPt+1)<< "\t"<< h_mean_p_betaReso[pid][iEta]->GetBinContent(iPt+1) << std::endl;
          DrawGammaLines(h_mean_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1), 
                        h_mean_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1), 
                        0., 0.05*h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kGray+2, 2);
          DrawGammaLines(h_mean_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1)-h_sigma_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1), 
                        h_mean_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1)-h_sigma_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1), 
                        0., 0.05*h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);
          DrawGammaLines(h_mean_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1)+h_sigma_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1), 
                        h_mean_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1)+h_sigma_p_betaSmearT0Reso[pid][iEta]->GetBinContent(iPt+1), 
                        0., 0.05*h_betaSmearT0_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);

          padnum++;
        }
        cPNG->Print(Form("%s/Beta_%s_ResolutionWithT0_%d_%d.%s", outputDirBetaRes.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      }
    }
  }
  cout << __LINE__ << endl;
  
  split_canvas(cPNG, "cPNG1", nEta+1);
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for(Int_t iEta = 0; iEta < nEta+1; iEta++){
      cPNG->cd(padnum+1);
      cPNG->cd(padnum+1)->Clear();
      padnum++;
    }
    padnum =0;
    
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.11, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      cPNG->cd(padnum+1)->SetLogx(1);
      SetStyleHistoTH2ForGraphs(h_beta_reso[pid][iEta], "#it{p} (GeV/#it{c})", "(1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(h_beta_reso[pid][iEta]);
      h_beta_reso[pid][iEta]->GetYaxis()->SetRangeUser(-0.02,0.02);
      h_beta_reso[pid][iEta]->GetZaxis()->SetRangeUser(0.9/5,h_beta_reso[pid][iEta]->GetMaximum()*2);
      h_beta_reso[pid][iEta]->GetXaxis()->SetRangeUser(0.1, 50);
      h_beta_reso[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
    }
    cPNG->Print(Form("%s/BetaP_%s_Resolution.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
  }
  
    
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for(Int_t iEta = 0; iEta < nEta+1; iEta++){
      cPNG->cd(padnum+1);
      cPNG->cd(padnum+1)->Clear();
      padnum++;
    }
    padnum =0;
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      if (!enablePIDEtaSmear[pid][iEta]) continue;
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.11, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      cPNG->cd(padnum+1)->SetLogx(1);
      SetStyleHistoTH2ForGraphs(h_betaSmearT0_reso[pid][iEta], "#it{p} (GeV/#it{c})", "(1/#beta^{rec}-1/#beta^{MC})/1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(h_betaSmearT0_reso[pid][iEta]);
      h_betaSmearT0_reso[pid][iEta]->GetYaxis()->SetRangeUser(-0.02,0.02);
      h_betaSmearT0_reso[pid][iEta]->GetXaxis()->SetRangeUser(0,50.);
      h_betaSmearT0_reso[pid][iEta]->GetZaxis()->SetRangeUser(0.9/5,h_betaSmearT0_reso[pid][iEta]->GetMaximum()*2);
      h_betaSmearT0_reso[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iEta == 0)  drawLatexAdd("w/ t_{0}",0.16,0.86-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
    }
    cPNG->Print(Form("%s/BetaP_%s_ResolutionWithT0.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
  }
  cout << __LINE__ << endl;
  

  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;

    for(Int_t iEta = 0; iEta < nEta+1; iEta++){
      cPNG->cd(padnum+1);
      cPNG->cd(padnum+1)->Clear();
      padnum++;
    }
    padnum =0;

    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      cPNG->cd(padnum+1)->SetLogx(1);
      SetStyleHistoTH2ForGraphs(h_betaGen_p[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(h_betaGen_p[pid][iEta]);
//       h_betaGen_p[pid][iEta]->GetYaxis()->SetRangeUser(-0.02,0.02);
      h_betaGen_p[pid][iEta]->GetZaxis()->SetRangeUser(0.9/5,h_betaGen_p[pid][iEta]->GetMaximum()*2);
      h_betaGen_p[pid][iEta]->GetXaxis()->SetRangeUser(0.1,50.);
      h_betaGen_p[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
      
      cSingle2D->cd();
        h_betaGen_p[pid][iEta]->Draw("colz");
  
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pid != 0) drawLatexAdd(Form("%s", partLabel[pid].Data()),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/GenBetaP_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/GenBetaP_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
  
  cout << __LINE__ << endl;
  
  for (Int_t pid = 0; pid < 6; pid++){
    for (Int_t eR = 0; eR < 3; eR++){
      if (!h_betaGen_p_Region[pid][eR]) continue;
      cSingle2D->cd();
      cSingle2D->SetLogx(1);
      SetStyleHistoTH2ForGraphs(h_betaGen_p_Region[pid][eR], "#it{p} (GeV/#it{c})", "1/#beta^{MC}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      
        ScaleByBinWidth2D(h_betaGen_p_Region[pid][eR]);
  //       h_betaGen_p_Region[pid][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
        h_betaGen_p_Region[pid][eR]->GetZaxis()->SetRangeUser(0.9/5,h_betaGen_p_Region[pid][eR]->GetMaximum()*2);
        h_betaGen_p_Region[pid][eR]->GetXaxis()->SetRangeUser(0.1,50.);
        h_betaGen_p_Region[pid][eR]->Draw("colz");
      
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]), 0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pid != 0) drawLatexAdd(Form("%s", partLabel[pid].Data()),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/GenBetaP_%s_Bin_%s.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data()));
    }
    cSingle2D->SetLogx(kFALSE);
  }
   
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for(Int_t iEta = 0; iEta < nEta+1; iEta++){
      cPNG->cd(padnum+1);
      cPNG->cd(padnum+1)->Clear();
      padnum++;
    }
    padnum =0;
    
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      cPNG->cd(padnum+1)->SetLogx(1);
      SetStyleHistoTH2ForGraphs(h_betaTrack_p[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{path}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(h_betaTrack_p[pid][iEta]);
      h_betaTrack_p[pid][iEta]->GetXaxis()->SetRangeUser(0,50.);
      h_betaTrack_p[pid][iEta]->GetZaxis()->SetRangeUser(0.9/5,h_betaTrack_p[pid][iEta]->GetMaximum()*2);
      h_betaTrack_p[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;

      cSingle2D->cd();
        h_betaTrack_p[pid][iEta]->Draw("colz");
        
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pid != 0) drawLatexAdd(Form("%s", partLabel[pid].Data()),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/TrackBetaP_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    }
    cPNG->Print(Form("%s/TrackBetaP_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
  cout << __LINE__ << endl;
  
  for (Int_t pid = 0; pid < 6; pid++){
    for (Int_t eR = 0; eR < 3; eR++){
      if (!h_betaTrack_p_Region[pid][eR]) continue;
      cSingle2D->cd();
      cSingle2D->SetLogx();
      SetStyleHistoTH2ForGraphs(h_betaTrack_p_Region[pid][eR], "#it{p} (GeV/#it{c})",  "1/#beta^{path}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(h_betaTrack_p_Region[pid][eR]);
  //       h_betaTrack_p_Region[pid][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
        h_betaTrack_p_Region[pid][eR]->GetZaxis()->SetRangeUser(0.9/5,h_betaTrack_p_Region[pid][eR]->GetMaximum()*2);
        h_betaTrack_p_Region[pid][eR]->GetXaxis()->SetRangeUser(0.07,50.);
        h_betaTrack_p_Region[pid][eR]->Draw("colz");
      
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pid != 0) drawLatexAdd(Form("%s", partLabel[pid].Data()),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/TrackBetaP_%s_Bin_%s.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data()));
    }
    cSingle2D->SetLogx(kFALSE);
  }
  
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for(Int_t iEta = 0; iEta < nEta+1; iEta++){
      cPNG->cd(padnum+1);
      cPNG->cd(padnum+1)->Clear();
      padnum++;
    }
    padnum =0;
    
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      cPNG->cd(padnum+1)->SetLogx(1);
      SetStyleHistoTH2ForGraphs(h_beta_p[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{rec}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(h_beta_p[pid][iEta]);
      h_beta_p[pid][iEta]->GetXaxis()->SetRangeUser(0.1,50);
      h_beta_p[pid][iEta]->GetZaxis()->SetRangeUser(0.9/5, h_beta_p[pid][iEta]->GetMaximum()*2);
      h_beta_p[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;

      cSingle2D->cd();
      cSingle2D->SetLogx(1);
        h_beta_p[pid][iEta]->Draw("colz");
        
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pid != 0) drawLatexAdd(Form("%s", partLabel[pid].Data()),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    

      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/BetaP_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/BetaP_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
  cout << __LINE__ << endl;
  for (Int_t pid = 0; pid < 6; pid++){
    for (Int_t eR = 0; eR < 3; eR++){
      if (!h_beta_p_Region[pid][eR]) continue;
      cSingle2D->cd();
      cSingle2D->SetLogx();
      SetStyleHistoTH2ForGraphs(h_beta_p_Region[pid][eR], "#it{p} (GeV/#it{c})",  "1/#beta^{rec}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
        ScaleByBinWidth2D(h_beta_p_Region[pid][eR]);
  //       h_beta_p_Region[pid][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
        h_beta_p_Region[pid][eR]->GetZaxis()->SetRangeUser(0.9/5,h_beta_p_Region[pid][eR]->GetMaximum()*2);
        h_beta_p_Region[pid][eR]->GetXaxis()->SetRangeUser(0.07,50.);
        h_beta_p_Region[pid][eR]->Draw("colz");
      
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pid != 0) drawLatexAdd(Form("%s", partLabel[pid].Data()),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/BetaP_%s_Bin_%s.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data()));
    }
    cSingle2D->SetLogx(kFALSE);
  }
  
  cout << __LINE__ << endl;
  for (Int_t pid = 0; pid < 6; pid++){
    Int_t padnum = 0;
    for(Int_t iEta = 0; iEta < nEta+1; iEta++){
      cPNG->cd(padnum+1);
      cPNG->cd(padnum+1)->Clear();
      padnum++;
    }
    padnum =0;
    
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      if (!enablePIDEtaSmear[pid][iEta]) continue;
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }      
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.095, 0.03, 0.045, 0.115);
      cPNG->cd(padnum+1)->SetLogz(1);
      cPNG->cd(padnum+1)->SetLogx(1);
      SetStyleHistoTH2ForGraphs(h_betaSmearT0_p[pid][iEta], "#it{p} (GeV/#it{c})", "1/#beta^{rec}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(h_betaSmearT0_p[pid][iEta]);
      h_betaSmearT0_p[pid][iEta]->GetXaxis()->SetRangeUser(0.1,50.);
      h_betaSmearT0_p[pid][iEta]->GetZaxis()->SetRangeUser(0.9/5,h_betaSmearT0_p[pid][iEta]->GetMaximum()*2);
      h_betaSmearT0_p[pid][iEta]->Draw("col");
      
      if(iEta == 0)  drawLatexAdd(Form("%s", partLabel[pid].Data()),0.16,0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iEta == 0)  drawLatexAdd("w/ t_{0}",0.16,0.86-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.85, 0.86,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      padnum++;
      
      cSingle2D->cd();
      cSingle2D->SetLogx(1);
        h_betaSmearT0_p[pid][iEta]->Draw("colz");
        
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pid != 0)
          drawLatexAdd(Form("w/ t_{0}, %s", partLabel[pid].Data()),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        else 
          drawLatexAdd("w/ t_{0}",0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/BetaPWithT0_%s_Bin_%d_%d.%s", outputDirBeta.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cPNG->Print(Form("%s/BetaPWithT0_%s.%s", outputDirBeta.Data(), partName[pid].Data(), suffix.Data()));
  }
 

  //*****************************************************************************
  // 1/beta plots in regions for different particle type
  //*****************************************************************************
  for (Int_t pid = 0; pid < 6; pid++){
    for (Int_t eR = 0; eR < 3; eR++){
      SetPlotStyle();
      if (!h_betaSmearT0_p_Region[pid][eR]) continue;
      cSingle2D->cd();
      cSingle2D->SetLogx();
      Double_t etaMin = partEta[minEtaBin[eR]];
      Double_t etaMax = partEta[maxEtaBin[eR]+1];
      SetStyleHistoTH2ForGraphs(h_betaSmearT0_p_Region[pid][eR], "#it{p} (GeV/#it{c})",  "1/#beta^{rec}", 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
        ScaleByBinWidth2D(h_betaSmearT0_p_Region[pid][eR]);
  //       h_betaSmearT0_p_Region[pid][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
        h_betaSmearT0_p_Region[pid][eR]->GetZaxis()->SetRangeUser(0.9/5,h_betaSmearT0_p_Region[pid][eR]->GetMaximum()*5);
        h_betaSmearT0_p_Region[pid][eR]->GetXaxis()->SetRangeUser(0.07,50.);
        h_betaSmearT0_p_Region[pid][eR]->GetYaxis()->SetRangeUser(0.9,1.5);
        h_betaSmearT0_p_Region[pid][eR]->Draw("colz");
     
        if (pid == 0){
          TLegend* legendPIDSummary   = GetAndSetLegend2(0.85 - 0.07*4, 0.91-(nLinesCol+3)*0.85*textSizeLabelsRel, 0.85, 0.91-(nLinesCol+4)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel, 4, "", 42, 0.5);
          
          for (Int_t id = 1; id < 6; id++){
            if (id == 2) continue;
            DrawGammaSetMarkerTGraphErr(   g_betaMeanSmearT0_p_Region[id][eR], 0, 0, colorPID[id], colorPID[id], 3, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaMeanSmearT0_p_Region[id][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaUpSmearT0_p_Region[id][eR], 0, 0, colorPID[id]+1, colorPID[id]+1, 2, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaSigmaUpSmearT0_p_Region[id][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaDownSmearT0_p_Region[id][eR], 0, 0, colorPID[id]+1, colorPID[id]+1, 2, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaSigmaDownSmearT0_p_Region[id][eR]->Draw("same, cx");
            legendPIDSummary->AddEntry(g_betaMeanSmearT0_p_Region[id][eR], partLabel[id], "l");
          }
          legendPIDSummary->Draw();
        } else {
          if (g_betaMeanSmearT0_p_Region[pid][eR]){
            DrawGammaSetMarkerTGraphErr(   g_betaMeanSmearT0_p_Region[pid][eR], 0, 0, colorPID[pid], colorPID[pid], 3, kFALSE, 0, kFALSE, lineStylePID[pid]);
            g_betaMeanSmearT0_p_Region[pid][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaUpSmearT0_p_Region[pid][eR], 0, 0, colorPID[pid]+1, colorPID[pid]+1, 2, kFALSE, 0, kFALSE, lineStylePID[pid]);
            g_betaSigmaUpSmearT0_p_Region[pid][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaDownSmearT0_p_Region[pid][eR], 0, 0, colorPID[pid]+1, colorPID[pid]+1, 2, kFALSE, 0, kFALSE, lineStylePID[pid]);
            g_betaSigmaDownSmearT0_p_Region[pid][eR]->Draw("same, cx");
          }
        }
        
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.85,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pid != 0) 
          drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps, %s", tRes, partLabel[pid].Data()),0.85,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        else
          drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.85, 0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        
      cSingle2D->cd();
      cSingle2D->SaveAs(Form("%s/BetaPWithT0_%s_Bin_%s.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data()));
      
      if (pid == 0){
        gStyle->SetPalette(kGreyScale);
        TColor::InvertPalette();      
        
        cSingle2D->cd();
        cSingle2D->SetLogx();
          h_betaSmearT0_p_Region[pid][eR]->Draw("colz");
      
          TLegend* legendPIDSummary   = GetAndSetLegend2(0.85 - 0.07*4, 0.91-(nLinesCol+3)*0.85*textSizeLabelsRel, 0.85, 0.91-(nLinesCol+4)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel, 4, "", 42, 0.5);
          
          for (Int_t id = 1; id < 6; id++){
            if (id == 2) continue;
            DrawGammaSetMarkerTGraphErr(   g_betaMeanSmearT0_p_Region[id][eR], 0, 0, colorPID[id], colorPID[id], 3, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaMeanSmearT0_p_Region[id][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaUpSmearT0_p_Region[id][eR], 0, 0, colorPIDdarker[id], colorPIDdarker[id], 2, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaSigmaUpSmearT0_p_Region[id][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaDownSmearT0_p_Region[id][eR], 0, 0, colorPIDdarker[id], colorPIDdarker[id], 2, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaSigmaDownSmearT0_p_Region[id][eR]->Draw("same, cx");
            legendPIDSummary->AddEntry(g_betaMeanSmearT0_p_Region[id][eR], partLabel[id], "l");
          }
          legendPIDSummary->Draw();
        
          drawLatexAdd(Form("%1.1f < #eta < %1.1f",partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.85,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.85, 0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          
        cSingle2D->cd();
        cSingle2D->SaveAs(Form("%s/BetaPWithT0_%s_Bin_%s_BW.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data()));

        cSingle2D->SetLogx(0);
          h_betaSmearT0_p_Region[pid][eR]->GetXaxis()->SetRangeUser(0.,5.);
          h_betaSmearT0_p_Region[pid][eR]->Draw("colz");
      
          for (Int_t id = 1; id < 6; id++){
            if (id == 2) continue;
            DrawGammaSetMarkerTGraphErr(   g_betaMeanSmearT0_p_Region[id][eR], 0, 0, colorPID[id], colorPID[id], 3, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaMeanSmearT0_p_Region[id][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaUpSmearT0_p_Region[id][eR], 0, 0, colorPIDdarker[id], colorPIDdarker[id], 2, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaSigmaUpSmearT0_p_Region[id][eR]->Draw("same, cx");
            DrawGammaSetMarkerTGraphErr(   g_betaSigmaDownSmearT0_p_Region[id][eR], 0, 0, colorPIDdarker[id], colorPIDdarker[id], 2, kFALSE, 0, kFALSE, lineStylePID[id]);
            g_betaSigmaDownSmearT0_p_Region[id][eR]->Draw("same, cx");
          }
          legendPIDSummary->Draw();
        
          drawLatexAdd(Form("%1.1f < #eta < %1.1f",partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.85,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.85, 0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        cSingle2D->SaveAs(Form("%s/BetaPWithT0_%s_Bin_%s_BW_linX.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data()));

        
      }
    }
  }

    //*****************************************************************************
  // NSigma plots in regions for different particle type
  //*****************************************************************************
  for (Int_t pid = 0; pid < 6; pid++){
    for (Int_t pidref = 1; pidref < 6; pidref++){
      for (Int_t eR = 0; eR < 3; eR++){
        SetPlotStyle();
        if (!h_NSigma_p_Region[pid][pidref][eR]) continue;
        cSingle2D->cd();
        cSingle2D->SetLogx();
        Double_t etaMin = partEta[minEtaBin[eR]];
        Double_t etaMax = partEta[maxEtaBin[eR]+1];
        SetStyleHistoTH2ForGraphs(h_NSigma_p_Region[pid][pidref][eR], "#it{p} (GeV/#it{c})",  Form("n #sigma_{%s}",partLabel[pidref].Data()), 
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
          ScaleByBinWidth2D(h_NSigma_p_Region[pid][pidref][eR]);
    //       h_NSigma_p_Region[pid][pidref][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
          h_NSigma_p_Region[pid][pidref][eR]->GetZaxis()->SetRangeUser(0.9/5,h_NSigma_p_Region[pid][pidref][eR]->GetMaximum()*5);
          h_NSigma_p_Region[pid][pidref][eR]->GetXaxis()->SetRangeUser(0.07,50.);
          h_NSigma_p_Region[pid][pidref][eR]->GetYaxis()->SetRangeUser(-15,15);
          h_NSigma_p_Region[pid][pidref][eR]->Draw("colz");
      
          drawLatexAdd(Form("%1.1f < #eta < %1.1f",partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.85,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pid != 0) 
            drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps, %s", tRes, partLabel[pid].Data()),0.85,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          else
            drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.85, 0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          
        cSingle2D->cd();
        cSingle2D->SaveAs(Form("%s/NSigma%sTOFBetaPWithT0_%s_Bin_%s.%s", outputDirOverview.Data(), partName[pidref].Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data()));
      
        
      }
    }
  }
  cSingle2D->SetLogx(kFALSE);
  
  //*****************************************************************************
  // 1/beta plots in nominal acceptance
  //*****************************************************************************
  for (Int_t eR = 0; eR < 3; eR++){
    if (!hBetaVsPDet[eR]) continue;
    cSingle2D->cd();
    cSingle2D->SetLogx();
    SetStyleHistoTH2ForGraphs(hBetaVsPDet[eR], "#it{p} (GeV/#it{c})",  "1/#beta^{rec}", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
      ScaleByBinWidth2D(hBetaVsPDet[eR]);
//       h_betaSmearT0_p_Region[pid][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
      hBetaVsPDet[eR]->GetZaxis()->SetRangeUser(0.9/5,hBetaVsPDet[eR]->GetMaximum()*2);
      hBetaVsPDet[eR]->GetXaxis()->SetRangeUser(0.07,50.);
      hBetaVsPDet[eR]->Draw("colz");
    
      drawLatexAdd(Form("%1.1f < #eta < %1.1f",nominalEtaRegion[eR][0],nominalEtaRegion[eR][1]),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    cSingle2D->cd();
    cSingle2D->SaveAs(Form("%s/BetaPWithT0_All_NomAcc_%s.%s", outputDirOverview.Data(),nameOutEtaRange[eR].Data(), suffix.Data()));
  }
  cSingle2D->SetLogx(kFALSE);

  //********************************************************************************
  // 1/beta_rec - 1/beta gen plots in nominal acceptance different particles
  //********************************************************************************
  DrawGammaCanvasSettings( cSingle2D, 0.115, 0.12, 0.025, 0.105);
  for (Int_t eR = 0; eR < 3; eR++){
    for (Int_t pid = 1; pid < 6; pid++){
      if (!hResBetaVsPDet_woT0[eR][pid]) continue;
      cSingle2D->cd();
      cSingle2D->SetLogx();
      SetStyleHistoTH2ForGraphs(hResBetaVsPDet_woT0[eR][pid], "#it{p} (GeV/#it{c})",  Form("1/#beta^{rec} - 1/#beta^{rec, %s exp.}", partLabel[pid].Data()),
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.05, 510, 505);
        ScaleByBinWidth2D(hResBetaVsPDet_woT0[eR][pid]);
  //       h_betaSmearT0_p_Region[pid][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
        hResBetaVsPDet_woT0[eR][pid]->GetZaxis()->SetRangeUser(0.9/5,hResBetaVsPDet_woT0[eR][pid]->GetMaximum()*2);
        hResBetaVsPDet_woT0[eR][pid]->GetXaxis()->SetRangeUser(0.07,50.);
        hResBetaVsPDet_woT0[eR][pid]->Draw("colz");
      
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",nominalEtaRegion[eR][0],nominalEtaRegion[eR][1]),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      cSingle2D->SaveAs(Form("%s/BetaResP_NomAcc_%s_%s.%s", outputDirOverview.Data(),nameOutEtaRange[eR].Data(), partName[pid].Data(), suffix.Data()));
    }
  }
    
  //********************************************************************************
  // 1/beta_rec with T0 - 1/beta gen plots in nominal acceptance different particles
  //********************************************************************************
  for (Int_t eR = 0; eR < 3; eR++){
    for (Int_t pid = 1; pid < 6; pid++){
      if (!hResBetaVsPDet_wT0[eR][pid]) continue;
      cSingle2D->cd();
      cSingle2D->SetLogx();
      SetStyleHistoTH2ForGraphs(hResBetaVsPDet_wT0[eR][pid], "#it{p} (GeV/#it{c})",  Form("1/#beta^{rec} - 1/#beta^{rec, %s exp.}", partLabel[pid].Data()), 
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.05, 510, 505);
        ScaleByBinWidth2D(hResBetaVsPDet_wT0[eR][pid]);
  //       h_betaSmearT0_p_Region[pid][eR]->GetYaxis()->SetRangeUser(-0.02,0.02);
        hResBetaVsPDet_wT0[eR][pid]->GetZaxis()->SetRangeUser(0.9/5,hResBetaVsPDet_wT0[eR][pid]->GetMaximum()*2);
        hResBetaVsPDet_wT0[eR][pid]->GetXaxis()->SetRangeUser(0.07,50.);
        hResBetaVsPDet_wT0[eR][pid]->Draw("colz");
      
        drawLatexAdd(Form("%1.1f < #eta < %1.1f",nominalEtaRegion[eR][0],nominalEtaRegion[eR][1]),0.14,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ), 0.14,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.85,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem,0.85,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.91-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      cSingle2D->SaveAs(Form("%s/BetaResPWithT0_NomAcc_%s_%s.%s", outputDirOverview.Data(),nameOutEtaRange[eR].Data(), partName[pid].Data(), suffix.Data()));
    }
  }


  //*****************************************************************************
  // Setup canvas & dummy hists for 1D Plots
  //*****************************************************************************
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
  cout << __LINE__ << endl;
  
  //*****************************************************************************
  // Plot PID resolution different eta slices for single particle
  //*****************************************************************************
  for (Int_t pid = 0; pid < 6; pid++){
    //**************************************************************************
    // mean single particle diff eta
    //**************************************************************************
    histoDummyBetaResMean->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_mean_p_betaReso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_mean_p_betaReso[pid][iEta]->Draw("same,p");
      if (pid == 0){
        if (iEta == nEta )
          legendPtResM->AddEntry(h_mean_p_betaReso[pid][iEta],Form("%1.1f < #eta < %1.1f",partEta[0],partEta[iEta]),"p");
        else 
          legendPtResM->AddEntry(h_mean_p_betaReso[pid][iEta],Form("%1.1f < #eta < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
    }
    legendPtResM->Draw();
    DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);
   
    drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s",partLabel[pid].Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolution_%s_Mean_p.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
    
    //**************************************************************************
    // sigma single particle diff eta
    //**************************************************************************
    histoDummyBetaResSigma->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_sigma_p_betaReso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_sigma_p_betaReso[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s",partLabel[pid].Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolution_%s_Sigma_p.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
    
    //**************************************************************************
    // mean single particle diff eta with T0
    //**************************************************************************
    histoDummyBetaResMean->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_mean_p_betaSmearT0Reso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_mean_p_betaSmearT0Reso[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);

    drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s w/ t_{0}",partLabel[pid].Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithT0_%s_Mean_p.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
    
    //**************************************************************************
    // sigma single particle diff eta regions with T0
    //**************************************************************************
    for (Int_t eR = 0; eR < 3; eR++){
      TLegend* legendPtResMeanRegion  = GetAndSetLegend2(0.14, 0.91-(maxNEtaBinsFull[eR]/2*0.85*textSizeLabelsRel), 0.6, 0.91,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      
      histoDummyBetaResMean->Draw();
      for(Int_t iEta=minEtaBinFull[eR]; iEta<maxEtaBinFull[eR];iEta++){
        if (!enablePIDEta[pid][iEta]) continue;
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_mean_p_betaSmearT0Reso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_mean_p_betaSmearT0Reso[pid][iEta]->Draw("same,p");
        legendPtResMeanRegion->AddEntry(h_mean_p_betaSmearT0Reso[pid][iEta],Form("%1.1f < #eta < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
      legendPtResMeanRegion->Draw();
      drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      if (pid != 0) 
        drawLatexAdd(Form("w/ t_{0}, %s",partLabel[pid].Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      else 
        drawLatexAdd("w/ t_{0}",0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/BetaResolutionWithT0_%s_Mean_p_%s.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data() ));
    }

    
    //**************************************************************************
    // sigma single particle diff eta with T0
    //**************************************************************************
    histoDummyBetaResSigma->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePIDEta[pid][iEta]) continue;
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_sigma_p_betaSmearT0Reso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_sigma_p_betaSmearT0Reso[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      if (pid != 0) 
        drawLatexAdd(Form("w/ t_{0}, %s",partLabel[pid].Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      else 
        drawLatexAdd("w/ t_{0}",0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithT0_%s_Sigma_p.%s", outputDirBetaRes.Data(), partName[pid].Data(), suffix.Data()));
    
    //**************************************************************************
    // sigma single particle diff eta regions with T0
    //**************************************************************************
    for (Int_t eR = 0; eR < 3; eR++){
      TLegend* legendPtResSigmaRegion  = GetAndSetLegend2(0.14, 0.91-(maxNEtaBinsFull[eR]/2*0.85*textSizeLabelsRel), 0.6, 0.91,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      
      histoDummyBetaResSigma->Draw();
      for(Int_t iEta=minEtaBinFull[eR]; iEta<maxEtaBinFull[eR];iEta++){
        if (!enablePIDEta[pid][iEta]) continue;
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_sigma_p_betaSmearT0Reso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_sigma_p_betaSmearT0Reso[pid][iEta]->Draw("same,p");
        legendPtResSigmaRegion->AddEntry(h_sigma_p_betaSmearT0Reso[pid][iEta],Form("%1.1f < #eta < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
      legendPtResSigmaRegion->Draw();
      drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      if (pid != 0) 
        drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps, %s",tRes, partLabel[pid].Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      else 
        drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/BetaResolutionWithT0_%s_Sigma_p_%s.%s", outputDirOverview.Data(), partName[pid].Data(), nameOutEtaRange[eR].Data(), suffix.Data() ));
    }
    
  }
  
  //*****************************************************************************
  // Plot PID resolution diff species together eta slices
  //*****************************************************************************
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
      if (!enablePIDEta[pid][iEta]) continue;
      DrawGammaSetMarker(h_sigma_p_betaReso[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_sigma_p_betaReso[pid][iEta]->Draw("same,p");      
      legendPtResPID->AddEntry(h_sigma_p_betaReso[pid][iEta],partLabel[pid].Data(),"p");
    }
    legendPtResPID->Draw();

    drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionPID_Sigma_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    histoDummyBetaResSigma->Draw();
    for (Int_t pid =0; pid < 6; pid++){
      if (!enablePIDEta[pid][iEta]) continue;
      DrawGammaSetMarker(h_sigma_p_betaSmearT0Reso[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_sigma_p_betaSmearT0Reso[pid][iEta]->Draw("same,p");      
    }
    legendPtResPID->Draw();
    drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.95,0.88-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithT0PID_Sigma_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    histoDummyBetaResMean->Draw();
    for (Int_t pid =0; pid < 6; pid++){
      if (!enablePIDEta[pid][iEta]) continue;
      DrawGammaSetMarker(h_mean_p_betaReso[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_mean_p_betaReso[pid][iEta]->Draw("same,p");      
    }
    legendPtResPID->Draw();

    drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionPID_Mean_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    histoDummyBetaResMean->Draw();
    for (Int_t pid =0; pid < 6; pid++){
      if (!enablePIDEta[pid][iEta]) continue;
      DrawGammaSetMarker(h_mean_p_betaSmearT0Reso[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_mean_p_betaSmearT0Reso[pid][iEta]->Draw("same,p");      
    }
    legendPtResPID->Draw();
    drawLatexAdd(perfLabel,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-2*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f < #eta < %1.1f",etaMin,etaMax),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("w/ t_{0}, #Deltat = %2.0f ps", tRes ),0.95,0.88-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/BetaResolutionWithT0PID_Mean_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  }

  //*****************************************************************************
  // plot T0 as function of nTracks
  //*****************************************************************************
  TCanvas* cTrackT0 = new TCanvas("cTrackT0","",0,0,1800,800);
  split_canvas(cTrackT0, "cTrackT0", 2, true);
  Int_t padnum = 0;
  Double_t maxNTracksZ = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    if (h_InitT0VsNtrk[tl]->GetMaximum() > maxNTracksZ) maxNTracksZ = h_InitT0VsNtrk[tl]->GetMaximum();
  }
  
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    DrawVirtualPadSettings( cTrackT0->cd(padnum+1), 0.11, 0.02, 0.03, 0.115);
    cTrackT0->cd(padnum+1)->SetLogz(1);
    SetStyleHistoTH2ForGraphs(h_InitT0VsNtrk[tl], "#it{N}_{trk}", "#it{t}_{0, initial} (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_InitT0VsNtrk[tl]->Scale(1./nEvents);
    h_InitT0VsNtrk[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_InitT0VsNtrk[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZ*1.2/nEvents);
    h_InitT0VsNtrk[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_InitT0VsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_InitT0VsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/InitialT0_NTrack.%s", outputDir.Data(), suffix.Data()));

  padnum = 0;
  Double_t maxNTracksZDiff = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    if (h_InitT0DiffVsNtrk[tl]->GetMaximum() > maxNTracksZDiff) maxNTracksZDiff = h_InitT0DiffVsNtrk[tl]->GetMaximum();
  }

  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_InitT0DiffVsNtrk[tl], "#it{N}_{trk}", "#it{t}_{0, initial} - #it{t}_{0, MC} (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_InitT0DiffVsNtrk[tl]->Scale(1./nEvents);
    h_InitT0DiffVsNtrk[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_InitT0DiffVsNtrk[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZDiff*1.2/nEvents);
    h_InitT0DiffVsNtrk[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 1)  drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.88-textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_InitT0DiffVsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_InitT0DiffVsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/InitialT0Diff_NTrack.%s", outputDir.Data(), suffix.Data()));

  padnum = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_T0VsNtrk[tl], "#it{N}_{trk}", "#it{t}_{0} (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_T0VsNtrk[tl]->Scale(1./nEvents);
    h_T0VsNtrk[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_T0VsNtrk[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZ*1.2/nEvents);
    h_T0VsNtrk[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 1)  drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.88-textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_T0VsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_T0VsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/T0_NTrack.%s", outputDir.Data(), suffix.Data()));

  padnum = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_T0DiffVsNtrk[tl], "#it{N}_{trk}", "#it{t}_{0} - #it{t}_{0, MC} (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_T0DiffVsNtrk[tl]->Scale(1./nEvents);
    h_T0DiffVsNtrk[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_T0DiffVsNtrk[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZDiff*1.2/nEvents);
    h_T0DiffVsNtrk[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 1)  drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.88-textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_T0DiffVsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_T0DiffVsNtrk[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/T0Diff_NTrack.%s", outputDir.Data(), suffix.Data()));
  
  //*****************************************************************************
  // plot T0 as function of nTracks with TTL
  //*****************************************************************************
  padnum = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_InitT0VsNttl[tl], "#it{N}_{trk w/ TTL}", "#it{t}_{0, initial} (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_InitT0VsNttl[tl]->Scale(1./nEvents);
    h_InitT0VsNttl[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_InitT0VsNttl[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZ*1.2/nEvents);
    h_InitT0VsNttl[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 1)  drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.88-textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_InitT0VsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_InitT0VsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/InitialT0_NTrackTTL.%s", outputDir.Data(), suffix.Data()));

  padnum = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_InitT0DiffVsNttl[tl], "#it{N}_{trk w/ TTL}", "#it{t}_{0, initial} - #it{t}_{0, MC}  (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_InitT0DiffVsNttl[tl]->Scale(1./nEvents);
    h_InitT0DiffVsNttl[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_InitT0DiffVsNttl[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZDiff*1.2/nEvents);
    h_InitT0DiffVsNttl[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 1)  drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.88-textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_InitT0DiffVsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_InitT0DiffVsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/InitialT0Diff_NTrackTTL.%s", outputDir.Data(), suffix.Data()));
  
  padnum = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_T0VsNttl[tl], "#it{N}_{trk w/ TTL}", "#it{t}_{0} (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_T0VsNttl[tl]->Scale(1./nEvents);
    h_T0VsNttl[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_T0VsNttl[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZ*1.2/nEvents);
    h_T0VsNttl[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 1)  drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.88-textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_T0VsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_T0VsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/T0_NTrackTTL.%s", outputDir.Data(), suffix.Data()));

  padnum = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h_T0DiffVsNttl[tl], "#it{N}_{trk w/ TTL}", "#it{t}_{0} - #it{t}_{0, MC} (ns)", 
                              0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
    h_T0DiffVsNttl[tl]->Scale(1./nEvents);
    h_T0DiffVsNttl[tl]->GetXaxis()->SetRangeUser(0,50.);
    h_T0DiffVsNttl[tl]->GetZaxis()->SetRangeUser(0.9/nEvents,maxNTracksZDiff*1.2/nEvents);
    h_T0DiffVsNttl[tl]->Draw("col");
    
    if(tl == 1)  drawLatexAdd(perfLabel,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 1)  drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.88-textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0)  drawLatexAdd(collisionSystem,0.93,0.88,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(tl == 0 && pTHard.CompareTo("") != 0 )  drawLatexAdd(pTHard,0.93,0.88-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
    if (tl == 0)drawLatexAdd(Form("%1.2f%s events w/ scat. e^{-}", h_T0DiffVsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if (tl == 1)drawLatexAdd(Form("%1.2f%s events w/o scat. e^{-}", h_T0DiffVsNttl[tl]->GetEntries()/nEvents*100,"%"),0.93, 0.16,textSizeSinglePad,kFALSE,kFALSE,kTRUE);

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/T0Diff_NTrackTTL.%s", outputDir.Data(), suffix.Data()));

  //*****************************************************************************
  // plot T0 initial different iterations
  //*****************************************************************************
  padnum = 0;
  for (Int_t tl = 0; tl < 2; tl++){
    cTrackT0->cd(padnum+1);
    DrawVirtualPadSettings( cTrackT0->cd(padnum+1), 0.1, 0.022, 0.02, 0.1);
    cTrackT0->cd(padnum+1)->SetLogz(0);
    cTrackT0->cd(padnum+1)->SetLogy(1);
    
    TLegend* legendIterT0 = GetAndSetLegend2(0.14, 0.905, 0.35, 0.905-nIter*0.85*textSizeSinglePad ,0.85*textSizeLabelsRel, 1, "", 42, 0.25);
    
    for (Int_t it = 0; it < nIter; it++){
      if (it == 0){
        SetStyleHistoTH1ForGraphs(h_T0VsIterProj[tl*nIter+it],"#it{t}_{0} - #it{t}_{MC} (ns)", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
        DrawGammaSetMarker(h_T0VsIterProj[tl*nIter+it], styleMarkerIter[it], sizeMarkerIter[it], colorIter[it], colorIter[it]);      
//         h_T0VsIterProj[tl]->GetXaxis()->SetRangeUser(0,50.);
        h_T0VsIterProj[tl*nIter+it]->GetYaxis()->SetRangeUser(10,nEvents/5);
        h_T0VsIterProj[tl*nIter+it]->Draw("hist");
        legendIterT0->AddEntry(h_T0VsIterProj[tl*nIter+it],Form("it %i, #sigma = %0.2f ps", it, h_T0VsIterProj[tl*nIter+it]->GetStdDev()*1000),"l");  
        
      } else if (it == nIter-1){
        DrawGammaSetMarker(h_T0VsIterProj[tl*nIter+it], styleMarkerIter[it], sizeMarkerIter[it], colorIter[it], colorIter[it]);      
        h_T0VsIterProj[tl*nIter+it]->Draw("pe,same");
        legendIterT0->AddEntry(h_T0VsIterProj[tl*nIter+it],Form("it %i, #sigma = %0.2f ps", it, h_T0VsIterProj[tl*nIter+it]->GetStdDev()*1000),"p");  
      } else {
        DrawGammaSetMarker(h_T0VsIterProj[tl*nIter+it], styleMarkerIter[it], sizeMarkerIter[it], colorIter[it], colorIter[it]);      
        h_T0VsIterProj[tl*nIter+it]->Draw("hist,same");
        legendIterT0->AddEntry(h_T0VsIterProj[tl*nIter+it],Form("it %i, #sigma = %0.2f ps", it, h_T0VsIterProj[tl*nIter+it]->GetStdDev()*1000),"l");  
      }
    }    
    h_T0VsIterProj[tl*nIter+0]->Draw("same,axis");
    if(tl == 1){
      drawLatexAdd(collisionSystem,0.93,0.91,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
      if(pTHard.CompareTo("") != 0 ){
        drawLatexAdd(pTHard,0.93,0.91-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(perfLabel,0.93,0.91-2*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
      } else {
         drawLatexAdd(perfLabel,0.93,0.91-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
         drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.93,0.91-2*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
      }
      drawLatexAdd("w/o scat. e^{-}",0.16, 0.91,textSizeSinglePad,kFALSE,kFALSE,kFALSE);
    }
    if (tl == 0){
      drawLatexAdd("w/ scat. e^{-}",0.16, 0.91,textSizeSinglePad,kFALSE,kFALSE,kFALSE);
    }
    legendIterT0->Draw();

    padnum++;
  }
  cTrackT0->SaveAs(Form("%s/InitialT0_Iterations.%s", outputDir.Data(), suffix.Data()));

  //*****************************************************************************
  // plot T0 as function of nTracks
  //*****************************************************************************
  TCanvas* cT0Ite = new TCanvas("cBetaSlice","",0,0,900,800);
  DrawGammaCanvasSettings( cT0Ite, 0.1, 0.022, 0.02, 0.1);
  cT0Ite->SetLogy(1);
  
  TLegend* legendIterT0 = GetAndSetLegend2(0.14, 0.905, 0.35, 0.905-nIter*0.85*textSizeSinglePad ,0.85*textSizeLabelsRel, 1, "", 42, 0.25);
    
    for (Int_t it = 0; it < nIter; it++){
      if (it == 0){
        SetStyleHistoTH1ForGraphs(h_T0VsIterProj[it],"#it{t}_{0} - #it{t}_{MC} (ns)", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
        DrawGammaSetMarker(h_T0VsIterProj[it], styleMarkerIter[it], sizeMarkerIter[it], colorIter[it], colorIter[it]);      
//         h_T0VsIterProj[tl]->GetXaxis()->SetRangeUser(0,50.);
        h_T0VsIterProj[it]->GetYaxis()->SetRangeUser(10,nEvents/5);
        h_T0VsIterProj[it]->Draw("hist");
        legendIterT0->AddEntry(h_T0VsIterProj[it],Form("it %i, #sigma = %0.2f ps", it, h_T0VsIterProj[it]->GetStdDev()*1000),"l");  
        
      } else if (it == nIter-1){
        DrawGammaSetMarker(h_T0VsIterProj[it], styleMarkerIter[it], sizeMarkerIter[it], colorIter[it], colorIter[it]);      
        h_T0VsIterProj[it]->Draw("pe,same");
        legendIterT0->AddEntry(h_T0VsIterProj[it],Form("it %i, #sigma = %0.2f ps", it, h_T0VsIterProj[it]->GetStdDev()*1000),"p");  
      } else {
        DrawGammaSetMarker(h_T0VsIterProj[it], styleMarkerIter[it], sizeMarkerIter[it], colorIter[it], colorIter[it]);      
        h_T0VsIterProj[it]->Draw("hist,same");
        legendIterT0->AddEntry(h_T0VsIterProj[it],Form("it %i, #sigma = %0.2f ps", it, h_T0VsIterProj[it]->GetStdDev()*1000),"l");  
      }
    }    
    h_T0VsIterProj[0]->Draw("same,axis");
    drawLatexAdd(collisionSystem,0.94,0.91,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    if(pTHard.CompareTo("") != 0 ){
      drawLatexAdd(pTHard,0.94,0.915-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(perfLabel,0.94,0.91-2*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    } else {
      drawLatexAdd(perfLabel,0.94,0.91-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#Deltat = %2.0f ps", tRes ),0.94,0.91-2*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
    }
    drawLatexAdd("w/ scat. e^{-}",0.16, 0.91,textSizeSinglePad,kFALSE,kFALSE,kFALSE);
    legendIterT0->Draw();
  cT0Ite->SaveAs(Form("%s/InitialT0_Iterations_WScatE.%s", outputDir.Data(), suffix.Data()));

  //*****************************************************************************
  // plot T0 as function of nTracks
  //*****************************************************************************
  TCanvas* cBetaSlice = new TCanvas("cBetaSlice","",0,0,900,800);
  DrawGammaCanvasSettings( cBetaSlice, 0.1, 0.022, 0.02, 0.1);
  cBetaSlice->SetLogy(1);
    
  for (Int_t eR = 0; eR <3 ; eR++){
    for (Int_t sl = 0; sl < nSlicesP; sl++){
      TLegend* legendSlice = GetAndSetLegend2(0.94 - 0.07*5, 0.91-(nLinesCol+4)*0.85*textSizeLabelsRel, 0.94, 0.91-(nLinesCol+5)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel, 5, "", 42, 0.5);
      for (Int_t pid = 0; pid < 6; pid++){
//         if (pid ==2 ) continue;
        hSliceBetaRegion[pid][eR][sl]->Rebin(4);
        if (pid == 0){
          SetStyleHistoTH1ForGraphs(hSliceBetaRegion[pid][eR][sl],"1/#it{#beta}", "counts",
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
          DrawGammaSetMarker(hSliceBetaRegion[pid][eR][sl], 0, 0, kGray+2, kGray+2);
          hSliceBetaRegion[pid][eR][sl]->SetLineWidth(2);
          hSliceBetaRegion[pid][eR][sl]->GetYaxis()->SetRangeUser(1,hSliceBetaRegion[pid][eR][sl]->GetMaximum()*2);
          hSliceBetaRegion[pid][eR][sl]->Draw("hist");
        } else {
          DrawGammaSetMarker(hSliceBetaRegion[pid][eR][sl], 0, 0, colorPID[pid], colorPID[pid]);
          hSliceBetaRegion[pid][eR][sl]->SetFillColorAlpha(colorPID[pid],0.5);
          hSliceBetaRegion[pid][eR][sl]->Draw("hist,same");
          legendSlice->AddEntry(hSliceBetaRegion[pid][eR][sl], partLabel[pid], "f");
        }
      }    
      hSliceBetaRegion[0][eR][sl]->Draw("same,hist");
      hSliceBetaRegion[0][eR][sl]->Draw("same,axis");
      
      drawLatexAdd(collisionSystem,0.93,0.91,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
      if(pTHard.CompareTo("") != 0 ){
        drawLatexAdd(pTHard,0.93,0.91-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(perfLabel,0.93,0.91-2*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(Form("%1.1f < #eta < %1.1f", partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.93,0.91-3*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(Form("p = %0.1f GeV", sliceP[sl]),0.93,0.91-4*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
      } else {
        drawLatexAdd(perfLabel,0.93,0.91-1*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(Form("%1.1f < #eta < %1.1f", partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]),0.93,0.91-2*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(Form("p = %0.1f GeV", sliceP[sl]),0.93,0.91-3*textSizeSinglePad,textSizeSinglePad,kFALSE,kFALSE,kTRUE);
      }
      legendSlice->Draw();  
      cBetaSlice->SaveAs(Form("%s/BetaExampleSlice_%s_%d.%s", outputDirOverview.Data(), nameOutEtaRange[eR].Data(), sl, suffix.Data()));
    }
  }
  
  textSizeLabelsPixel       = 55;
  textSizeLabelsRel      = 55./1200;

  TCanvas* canvasEoverPRejection       = new TCanvas("canvasEoverPRejection", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasEoverPRejection,  0.1, 0.01, 0.017, 0.095);
  canvasEoverPRejection->SetLogy(1);
  canvasEoverPRejection->SetLogx(1);

  TH2F * histo2DRejectionDummy;
  histo2DRejectionDummy                = new TH2F("histo2DRejectionDummy", "histo2DRejectionDummy",1000, 0.061, 3.1, 100, 0.9, 999999 );
  SetStyleHistoTH2ForGraphs( histo2DRejectionDummy, "#it{p}_{track} (GeV/#it{c})", "#pi^{-}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
  histo2DRejectionDummy->GetYaxis()->SetLabelOffset(0.001);
  histo2DRejectionDummy->GetXaxis()->SetNoExponent();
  histo2DRejectionDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DRejectionDummy->DrawCopy();
  
  Style_t markerStyleCalo[3] = {20, 47, 33};
  Size_t markerSizeCalo[3] = {1.5, 1.8, 2.2};
  Color_t colorCalo[3];
  colorCalo[0] = getCaloColor("EEMC", 0);
  colorCalo[1] = getCaloColor("BEMC", 0);
  colorCalo[2] = getCaloColor("FEMC", 0);
  TLegend* legendTTL   = GetAndSetLegend2(0.65, 0.95-(3*textSizeLabelsRel), 0.95, 0.95,0.9*textSizeLabelsPixel,1);
  
  for (Int_t eR = 0; eR < 3; eR++){
    shrinkGraphErr( g_ParticleRejection_p_Region[3][1][eR], 20);
    DrawGammaSetMarkerTGraphErr(g_ParticleRejection_p_Region[3][1][eR], markerStyleCalo[eR], 1.5*markerSizeCalo[eR], colorCalo[eR] , colorCalo[eR]);
    g_ParticleRejection_p_Region[3][1][eR]->Draw("samepe");
    legendTTL->AddEntry(g_ParticleRejection_p_Region[3][1][eR], Form("%1.1f < #eta < %1.1f", partEta[minEtaBinTTL[eR]],partEta[maxEtaBinTTL[eR]+1]), "p");
  }
  legendTTL->Draw();
  float toplegval = 0.90;
//   drawLatexAdd("#it{p}_{rec}",0.75,0.30,textSizeLabelsRel,false,false,false);
//   drawLatexAdd("#it{p}_{true}",0.65,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd(perfLabel.Data(),0.14,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.14,toplegval-0.05,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("AC-LGAD, #Deltat = %2.0f ps", tRes ),0.14,toplegval-0.05*2,textSizeLabelsRel,false,false,false);
  
  canvasEoverPRejection->Update();
  canvasEoverPRejection->Print(Form("%s/pionrejection_TOF.%s",outputDir.Data(),suffix.Data()));
  
  
  //*****************************************************************************
  // Write output
  //*****************************************************************************
  TFile* outputFile  = new TFile(inputFileName.Data(),"UPDATE");
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      if (properFit){
        for (Int_t pid = 0; pid < 6; pid++){
          if (!enablePIDEta[pid][iEta]) continue;
          h_mean_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_FitMean_%d", readTrackClass.Data(), partName[pid].Data(),iEta),TObject::kOverwrite);
          h_sigma_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_FitSigma_%d",readTrackClass.Data(),  partName[pid].Data(), iEta),TObject::kOverwrite);
          if (h_mean_p_betaSmearT0Reso[pid][iEta]) h_mean_p_betaSmearT0Reso[pid][iEta]->Write(Form("histBetaResol%sWT0%s_FitMean_%d", readTrackClass.Data(), partName[pid].Data(),iEta),TObject::kOverwrite);
          if (h_sigma_p_betaSmearT0Reso[pid][iEta]) h_sigma_p_betaSmearT0Reso[pid][iEta]->Write(Form("histBetaResol%sWT0%s_FitSigma_%d", readTrackClass.Data(),  partName[pid].Data(), iEta),TObject::kOverwrite);
        }
      } else {
        for (Int_t pid = 0; pid < 6; pid++){
          if (!enablePIDEta[pid][iEta]) continue;
          h_mean_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
          h_sigma_p_betaReso[pid][iEta]->Write(Form("histBetaResol%s_%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);      
          if (h_mean_p_betaSmearT0Reso[pid][iEta])h_mean_p_betaSmearT0Reso[pid][iEta]->Write(Form("histBetaResol%sWT0%s_mean_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
          if (h_sigma_p_betaSmearT0Reso[pid][iEta])h_sigma_p_betaSmearT0Reso[pid][iEta]->Write(Form("histBetaResol%sWT0%s_sigma_%d", readTrackClass.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);      
        }
      }
    }
    g_ParticleRejection_p_Region[3][1][0]->Write("PionRejection_ETTL",TObject::kOverwrite);
    g_ParticleRejection_p_Region[3][1][1]->Write("PionRejection_CTTL",TObject::kOverwrite);
    g_ParticleRejection_p_Region[3][1][2]->Write("PionRejection_FTTL",TObject::kOverwrite);
  outputFile->Write();
  outputFile->Close();

  TFile* outputFile2  = new TFile(Form("%s/TOFparams.root",outputDir.Data()),"UPDATE");
    for(Int_t eR = 0; eR < 3; eR++){
      g_ParticleRejection_p_Region[3][1][0]->Write(Form("PionRejection_%s",readRegion[eR].Data() ),TObject::kOverwrite);
      for (Int_t pid = 1; pid < 6; pid++){
        g_betaMeanSmearT0_p_Region[pid][eR]->Write(Form("meanBeta_%s_%s", partName[pid].Data(), readRegion[eR].Data()) ,TObject::kOverwrite);
        g_betaSigmaAbsSmearT0_p_Region[pid][eR]->Write(Form("simgaBeta_%s_%s", partName[pid].Data(), readRegion[eR].Data() ),TObject::kOverwrite);
        g_betaSigmaDownSmearT0_p_Region[pid][eR]->Write(Form("3simgaBetaDown_%s_%s", partName[pid].Data(), readRegion[eR].Data()) ,TObject::kOverwrite);
        g_betaSigmaUpSmearT0_p_Region[pid][eR]->Write(Form("3simgaBetaUp_%s_%s", partName[pid].Data(), readRegion[eR].Data()) ,TObject::kOverwrite);
      }
    }
  outputFile2->Write();
  outputFile2->Close();
  
}
