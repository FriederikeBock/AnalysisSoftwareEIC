#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void positionResolutionCalorimeters(
    TString inputFileName     = "output_ERH.root",
    TString caloName          = "BECAL",
    TString clusterizerName   = "MA",
    TString suffix            = "pdf",
    Bool_t doFitting          = kTRUE
){  

  
  
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem           = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString perfLabel                 = "ECCE G4 simulation";
  
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;

  Double_t minEtaMax  = -4;
  Double_t maxEtaMax  = 4;
  Double_t minPhiMax  = -4;
  Double_t maxPhiMax  = 4;
  Double_t maxReso    = 1.25;
  TString detLabel    = "";
  Bool_t isEMCal      = 0;
  Int_t rebinRes      = 1;
  Float_t maxResoPer  = 31;
  Float_t maxFitRes1oE  = -1;
  Float_t minFitRes1oE  = -1;
  Float_t maxFitResE    = -1;
  Float_t minFitResE    = -1;
  Float_t maxFitMeanE    = -1;
  Float_t minFitMeanE    = 10;
  TString caloNameRead = caloName;
  Bool_t isComb         = 0;
  
  Int_t exEbin          = 15;
  Int_t exEtabin        = 2;
  Int_t nbins           = 100; 

  
  // fwd Hcal numbers
  Float_t sqrtETerm[2]  = {35, 50};
  Float_t linTerm[2]    = {7, 10};
  
  
  
  if (caloName.CompareTo("EEMC") == 0){
    detLabel  = "EEMC (PbWO_{4} crystal)";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    isEMCal   = 1;
    maxResoPer  = 12.5;
    maxFitRes1oE  = 1.45;
    minFitRes1oE  = 1/TMath::Sqrt(20.);
    maxFitResE    = 20;
    minFitResE    = 0.2;
    exEtabin      = 2;
    sqrtETerm[0]  = 2;
    sqrtETerm[1]  = 3;
    linTerm[0]    = 1;
    linTerm[1]    = 3;
  } else if (caloName.CompareTo("FEMC") == 0){
    detLabel  = "FEMC (PHENIX re-use)";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    isEMCal   = 1;
    maxResoPer  = 20.5;
    exEtabin      = 11;
    sqrtETerm[0]  = 7;
    sqrtETerm[1]  = 10;
    linTerm[0]    = 1;
    linTerm[1]    = 3;
    maxFitMeanE   = 20;
  } else if (caloName.CompareTo("BECAL") == 0){
    detLabel  = "BECAL (Sci-glass)";
    minEtaMax = -1.8;
    maxEtaMax = 1.3;
    isEMCal   = 1;
    maxResoPer  = 17.5;

    exEtabin      = 7;
    sqrtETerm[0]  = 7;
    sqrtETerm[1]  = 10;
    linTerm[0]    = 1;
    linTerm[1]    = 3;
  } else if (caloName.CompareTo("LFHCAL") == 0){
    detLabel  = "LFHCAL";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    maxReso   = 1.65;
    rebinRes    = 2;
    maxResoPer  = 72;
    exEtabin      = 11;
  } else if (caloName.CompareTo("LFHCAL-wMat") == 0){
    detLabel  = "LFHCAL, FEMC infront";
    caloNameRead = "LFHCAL";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    maxReso   = 1.65;
    rebinRes    = 2;
    maxResoPer  = 72;
    exEtabin      = 11;
  } else if (caloName.CompareTo("LFHCAL-FEMC") == 0){
    detLabel  = "LFHCAL+FEMC";
    caloNameRead = "LFHCAL";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    maxReso   = 1.65;
    rebinRes    = 2;
    maxResoPer  = 72;
    isComb      = 1;
    exEtabin      = 11;
  } else if (caloName.CompareTo("EHCAL") == 0){
    detLabel  = "EHCAL (STAR re-use)";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    maxReso   = 1.25;
    maxResoPer  = 62;
    exEtabin      = 2;
    sqrtETerm[0]  = 45;
    sqrtETerm[1]  = 50;
    linTerm[0]    = 6;
    linTerm[1]    = 10;

  } else if (caloName.CompareTo("EHCAL-wMat") == 0){
    detLabel  = "EHCAL, EEMC infront";
    caloNameRead = "EHCAL";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    maxReso   = 1.25;
    maxResoPer  = 62;
    exEtabin      = 2;
    sqrtETerm[0]  = 45;
    sqrtETerm[1]  = 50;
    linTerm[0]    = 6;
    linTerm[1]    = 10;
  } else if (caloName.CompareTo("EHCAL-EEMC") == 0){
    detLabel  = "EHCAL+EEMC";
    caloNameRead = "EHCAL";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    maxReso   = 1.25;
    maxResoPer  = 62;
    isComb      = 1;
    exEtabin      = 2;
    sqrtETerm[0]  = 45;
    sqrtETerm[1]  = 50;
    linTerm[0]    = 6;
    linTerm[1]    = 10;
  } else if (caloName.CompareTo("CHCAL") == 0){
    detLabel  = "oHCAL, iHCal infront";
    caloNameRead = "HCALOUT";
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    maxReso   = 1.25;
    maxResoPer  = 82;
    isComb      = 0;
    exEtabin      = 2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
  } else if (caloName.CompareTo("CHCAL-comb") == 0){
    detLabel  = "oHCAL+iHCal";
    caloNameRead = "HCALOUT";
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    maxReso   = 1.25;
    maxResoPer  = 82;
    isComb      = 0;
    exEtabin      = 2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
    isComb      = 1;
  } else if (caloName.CompareTo("HCALIN") == 0){
    detLabel  = "iHCal";
    caloNameRead = "HCALIN";
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    maxReso   = 1.25;
    maxResoPer  = 82;
    isComb      = 0;
    exEtabin      = 2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
  } else {
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
  }
  // determine eta bins
  Int_t binEtaMin = 0;
  Int_t binEtaMax = 15;
  
  TString outputDir                 = Form("plotsResolution%s",caloName.Data());
  if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);


  //************************** Read data **************************************************
  const int nEne = 34;
  const int nEta = 40;
  const int nPhi = 25;
  
  const static Double_t partE[]   = {   0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5,
                                        10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 50, 67.5, 75, 100, 150
                                    };


  const static Double_t partPhi[]   = {   -3, 2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0,
                                          0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3
                                    };

  const static Double_t partEta[]   = {   -5, -4.75, -4.5, -4.25, -4, -3.75,  -3.5, -3.25, -3,
                                         -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0,
                                          0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 
                                          4.5, 4.75, 5
                                    };





  TGraphErrors* graphReq          = new TGraphErrors(nEne); 
  TGraphErrors* graphReq1oE       = new TGraphErrors(nEne); 

  graphReq1oE->Sort();
  DrawGammaSetMarkerTGraphErr(  graphReq, 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
  DrawGammaSetMarkerTGraphErr(  graphReq1oE, 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
  TLegend* legendResYR    = GetAndSetLegend2(0.14, 0.14, 0.6, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYR->AddEntry(graphReq, "YR requirement","f");
  legendResYR->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETerm[0], sqrtETerm[1], linTerm[0], linTerm[1]),"");
  TLegend* legendResYR1oE = GetAndSetLegend2(0.58, 0.14, 0.97, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYR1oE->AddEntry(graphReq1oE, "YR requirement","f");
  legendResYR1oE->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETerm[0], sqrtETerm[1], linTerm[0], linTerm[1]),"");
  
  
  Int_t useEta[7][nEta]             = {{ 0}};
  Int_t usePhi[7][nPhi]             = {{ 0}};
  Int_t usePart[7]                  = { 0};
  Int_t usePartPhi[7]               = { 0};
  Int_t usePartEta[7]               = { 0};

  Int_t useParEnetPhi[7]            = { 0};
  Int_t useEnePhi[7][nEne]         = {{ 0}};

  Int_t useParEnetEta[7]            = { 0};
  Int_t useEneEta[7][nEne]         = {{ 0}};
          

  Double_t fitrange_low[nEne]     = { 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                      0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6,0.65, 0.7, 
                                      0.7, 0.7, 0.7, 0.7};
  Double_t fitrange_hig[nEne]     = { 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 
                                      0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 
                                      0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 
                                      0.9, 0.9, 1.0, 1.0};

  Int_t nParticles                  = 7;
  TString readNames[7]              = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_", ""};
  TString labelPart[7]              = {"#gamma", "e^{#pm}", "#pi^{#pm}", "p/#bar{p}", "K^{#pm}", "n", "all"};
  Color_t colorPart[7]              = {807, kBlue+1, kRed+2, kGreen+2, kViolet+2, kGray+2, kBlack};
  Style_t markerStylePart[7]        = {20, 46, 33, 30, 28, 42, 25};
  Style_t lineStylePart[7]          = {3, 4, 5, 7, 9, 3, 1};
  Size_t markerSizePart[7]          = {2.5*1.5, 2*1.5, 3*1.5, 2.7*1.5, 2*1.5, 3*1.5, 2*1.5};
  
  Float_t minEPart[7];
  Float_t maxEPart[7];
  Int_t nActivePart                 = 0;
  
  Int_t nActiveEtapart[7]           = {0};

  Int_t nActivePhipart[7]           = {0};
  
  TFile* inputFile                      = new TFile(inputFileName.Data());
  TH2F* hist2DResoEta[7]                = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1F* hist1DResoEtabins[7][nEta]      = {{NULL}};
  TF1* fitResosEta[7][nEta]             = {{NULL}};
  TH1F* histResoEtaVsEta[7]             = {NULL};
  TH1F* histMeanEtaVsEta[7]             = {NULL};
  TH1F* histSigmaEtaVsEta[7]            = {NULL};

  TH1F* histResoPhiVsPhi[7]             = {NULL};
  TH1F* histMeanPhiVsPhi[7]             = {NULL};
  TH1F* histSigmaPhiVsPhi[7]            = {NULL};

  TH2F* hist2DResoPhi[7]                = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1F* hist1DResoPhibins[7][nEta]      = {{NULL}};
  TF1* fitResosPhi[7][nPhi]                    = {{NULL}};

  TH3F* hist3DResoEta[7]                      = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1F* hist1DResoEtabins3D[7][nEne][nEta]    = {{NULL}};
  TF1*  fitResosEEta[7][nEne][nEta]           = {{NULL}};
  TH1F* histResoEEtaVsEta[7][nEne]            = {{NULL}};
  TH1F* histMeanEEtaVsEta[7][nEne]            = {{NULL}};
  TH1F* histSigmaEEtaVsEta[7][nEne]           = {{NULL}};


  TH3F* hist3DResoPhi[7]                      = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1F* hist1DResoPhibins3D[7][nEne][nPhi]    = {{NULL}};
  TF1*  fitResosEPhi[7][nEne][nEta]           = {{NULL}};
  TH1F* histResoEPhiVsPhi[7][nEne]            = {{NULL}};
  TH1F* histMeanEPhiVsPhi[7][nEne]            = {{NULL}};
  TH1F* histSigmaEPhiVsPhi[7][nEne]           = {{NULL}};


 for (Int_t i = 0; i < nParticles; i++){

    Int_t nbins = 100; 
    Double_t mineta = -5.; Double_t maxeta = 5.; 

    TString tempName = Form("%s/h_ERH_Reso_%sE_Eta_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempName = Form("%s/h_RH_ResoComb_%sE_Eta_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );

    hist3DResoEta[i]    = (TH3F*)inputFile->Get(tempName.Data());

    if (hist3DResoEta[i]->GetEntries() > 10){

      for (int ebin = 0; ebin < nEne; ebin++){

        histResoEEtaVsEta[i][ebin]  = new TH1F(Form("h_%sreso_EEta_%i", readNames[i].Data(),ebin), "", nbins, mineta, maxeta);
        histMeanEEtaVsEta[i][ebin]  = new TH1F(Form("h_%sreso_EEtaMean_%i", readNames[i].Data(),ebin), "", nbins, mineta, maxeta);
        histSigmaEEtaVsEta[i][ebin] = new TH1F(Form("h_%sreso_EEtaSigma_%i", readNames[i].Data(),ebin), "", nbins, mineta, maxeta);


        for (int ebin2 = 0; ebin2 < nEta; ebin2++){

          hist1DResoEtabins3D[i][ebin][ebin2]    = (TH1F*)hist3DResoEta[i]->ProjectionZ(Form("projection_Eta_3D_%s%d", readNames[i].Data(), ebin), 
                                                                              hist3DResoEta[i]->GetXaxis()->FindBin(partE[ebin]-0.3), hist3DResoEta[i]->GetXaxis()->FindBin(partE[ebin]+0.3),
                                                                              hist3DResoEta[i]->GetYaxis()->FindBin(partEta[ebin2]-0.125),hist3DResoEta[i]->GetYaxis()->FindBin(partEta[ebin2]+0.125),"e");

          if (hist1DResoEtabins3D[i][ebin][ebin2]->GetEntries() > 10 ){

            useParEnetPhi[i]           = 1;
            useEnePhi[i][ebin]          = 1;

          
            Double_t mean     = FindLargestBin1DHist(hist1DResoEtabins3D[i][ebin][ebin2], -2, 2);
            Double_t meanErr  = hist1DResoEtabins3D[i][ebin][ebin2]->GetMeanError();
            Double_t sigma    = hist1DResoEtabins3D[i][ebin][ebin2]->GetRMS();
            Double_t sigmaErr = hist1DResoEtabins3D[i][ebin][ebin2]->GetRMSError();

            if (doFitting){    
              fitResosEEta[i][ebin][ebin2]     = new TF1(Form("fit_resoeeta_%s%d", readNames[i].Data(), ebin), "gaus", mean-2*sigma, mean+2*sigma);
              fitResosEEta[i][ebin][ebin2]->SetParLimits(1, 0.6*mean, 1.2*mean);
              fitResosEEta[i][ebin][ebin2]->SetParameter(1, mean);
              fitResosEEta[i][ebin][ebin2]->SetParLimits(2, 0.6*sigma,sigma);
              fitResosEEta[i][ebin][ebin2]->SetParameter(2, sigma);
              hist1DResoEtabins3D[i][ebin][ebin2]->Fit(fitResosEEta[i][ebin][ebin2],"L0RMEQ","",mean-1.5*sigma, mean+1.5*sigma);

              mean     = fitResosEEta[i][ebin][ebin2]->GetParameter(1);
              meanErr  = fitResosEEta[i][ebin][ebin2]->GetParError(1);
              sigma    = fitResosEEta[i][ebin][ebin2]->GetParameter(2);
              sigmaErr = fitResosEEta[i][ebin][ebin2]->GetParError(2);
            }
          
            Double_t reso     = sigma/mean;
            Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
            histResoEEtaVsEta[i][ebin]->SetBinContent(histResoEEtaVsEta[i][ebin]->FindBin(partEta[ebin2]),  100*reso);
            histResoEEtaVsEta[i][ebin]->SetBinError(histResoEEtaVsEta[i][ebin]->FindBin(partEta[ebin2]), resoErr);

            histMeanEEtaVsEta[i][ebin]->SetBinContent(histResoEEtaVsEta[i][ebin]->FindBin(partEta[ebin2]), mean);
            histMeanEEtaVsEta[i][ebin]->SetBinError(histResoEEtaVsEta[i][ebin]->FindBin(partEta[ebin2]), meanErr);
            histSigmaEEtaVsEta[i][ebin]->SetBinContent(histResoEEtaVsEta[i][ebin]->FindBin(partEta[ebin2]), sigma);
            histSigmaEEtaVsEta[i][ebin]->SetBinError(histResoEEtaVsEta[i][ebin]->FindBin(partEta[ebin2]), sigmaErr);  
          }  
          
        }
      }
    }
  }


  for (Int_t i = 0; i < nParticles-1; i++){

    if (!useParEnetPhi[i]) continue;

     TCanvas* phidistE = new TCanvas("phidistE","",0,0,1100,1000);
     DrawGammaCanvasSettings( phidistE, 0.095, 0.01, 0.01, 0.105);
     phidistE->Divide(1,2);
     phidistE->cd(1);
     
     Int_t pC = 0; 
     TLegend* legendEEta  = GetAndSetLegend2(0.2,0.3, 0.35, 0.95,0.4*textSizeLabelsPixel, 1, "", 43, 0.2);
     for (int ebin = 0; ebin < nEne; ebin++){
        if (!useEnePhi[i][ebin]) continue;
        
        drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsPixel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        
        if(useEnePhi[i][ebin]){

          SetStyleHistoTH1ForGraphs(histMeanEEtaVsEta[i][ebin], "#eta ^{MC}","(#eta ^{rec} - #eta ^{MC})_{Mean}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
          histMeanEEtaVsEta[i][ebin]->SetLineColor(pC+1);
          histMeanEEtaVsEta[i][ebin]->GetYaxis()->SetDecimals(kTRUE);
          histMeanEEtaVsEta[i][ebin]->GetYaxis()->SetMaxDigits(5);
          if(pC==0) {histMeanEEtaVsEta[i][ebin]->Draw("");}
          else  {histMeanEEtaVsEta[i][ebin]->Draw("SAME");}
          pC++;

          legendEEta->AddEntry(histMeanEEtaVsEta[i][ebin], Form("E^{MC} =  %1.1f", partE[ebin]),"lep");
        }   
     }
     legendEEta->Draw();

     phidistE->cd(2);
     
     pC = 0; 
     TLegend* legendEEtaSigma  = GetAndSetLegend2(0.2,0.3, 0.35, 0.95,0.4*textSizeLabelsPixel, 1, "", 43, 0.2);
     for (int ebin = 0; ebin < nEne; ebin++){
        if (!useEnePhi[i][ebin]) continue;
        if(useEnePhi[i][ebin]){
          SetStyleHistoTH1ForGraphs(histSigmaEEtaVsEta[i][ebin], "#eta ^{MC}","(#eta ^{rec} - #eta ^{MC})_{width}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
          histSigmaEEtaVsEta[i][ebin]->SetLineColor(pC+1);
          histSigmaEEtaVsEta[i][ebin]->GetYaxis()->SetDecimals(kTRUE);
          histSigmaEEtaVsEta[i][ebin]->GetYaxis()->SetMaxDigits(5);
          if(pC==0) {histSigmaEEtaVsEta[i][ebin]->Draw("");}
          else  {histSigmaEEtaVsEta[i][ebin]->Draw("SAME");}
          pC++;

          legendEEtaSigma->AddEntry(histSigmaEEtaVsEta[i][ebin], Form("E^{MC} =  %1.1f", partE[ebin]),"lep");
         }   
      }
      legendEEtaSigma->Draw();
      phidistE->Print(Form("%s/%setavsE_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));
  }


  for (Int_t i = 0; i < nParticles; i++){

    
    Double_t minphi = -5.; Double_t maxphi = 5.; 

    TString tempName = Form("%s/h_ERH_Reso_%sE_Phi_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempName = Form("%s/h_RH_ResoComb_%sE_Phi_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );

    hist3DResoPhi[i]    = (TH3F*)inputFile->Get(tempName.Data());

    if (hist3DResoPhi[i]->GetEntries() > 10){

      for (int ebin = 0; ebin < nEne; ebin++){

        histResoEPhiVsPhi[i][ebin]  = new TH1F(Form("h_%sreso_EPhi_%i", readNames[i].Data(),ebin), "", nbins, minphi, maxphi);
        histMeanEPhiVsPhi[i][ebin]  = new TH1F(Form("h_%sreso_EPhiMean_%i", readNames[i].Data(),ebin), "", nbins, minphi, maxphi);
        histSigmaEPhiVsPhi[i][ebin] = new TH1F(Form("h_%sreso_EPhiSigma_%i", readNames[i].Data(),ebin), "", nbins, minphi, maxphi);


        for (int ebin2 = 0; ebin2 < nPhi; ebin2++){

          hist1DResoPhibins3D[i][ebin][ebin2]    = (TH1F*)hist3DResoPhi[i]->ProjectionZ(Form("projection_Phi_3D_%s%d", readNames[i].Data(), ebin), 
                                                                              hist3DResoPhi[i]->GetXaxis()->FindBin(partE[ebin]-0.3), hist3DResoPhi[i]->GetXaxis()->FindBin(partE[ebin]+0.3),
                                                                              hist3DResoPhi[i]->GetYaxis()->FindBin(partPhi[ebin2]-0.125),hist3DResoPhi[i]->GetYaxis()->FindBin(partPhi[ebin2]+0.125),"e");

          if (hist1DResoPhibins3D[i][ebin][ebin2]->GetEntries() > 10 ){

            useParEnetPhi[i]        = 1;
            useEnePhi[i][ebin]      = 1;
          
            Double_t mean     = FindLargestBin1DHist(hist1DResoPhibins3D[i][ebin][ebin2], -2, 2);
            Double_t meanErr  = hist1DResoPhibins3D[i][ebin][ebin2]->GetMeanError();
            Double_t sigma    = hist1DResoPhibins3D[i][ebin][ebin2]->GetRMS();
            Double_t sigmaErr = hist1DResoPhibins3D[i][ebin][ebin2]->GetRMSError();

            if (doFitting){
                
              fitResosEPhi[i][ebin][ebin2]     = new TF1(Form("fit_resoeeta_%s%d", readNames[i].Data(), ebin), "gaus", mean-2*sigma, mean+2*sigma);
              fitResosEPhi[i][ebin][ebin2]->SetParLimits(1, 0.6*mean, 1.2*mean);
              fitResosEPhi[i][ebin][ebin2]->SetParameter(1, mean);
              fitResosEPhi[i][ebin][ebin2]->SetParLimits(2, 0.6*sigma,sigma);
              fitResosEPhi[i][ebin][ebin2]->SetParameter(2, sigma);
              hist1DResoPhibins3D[i][ebin][ebin2]->Fit(fitResosEPhi[i][ebin][ebin2],"L0RMEQ","",mean-10*sigma, mean+1.*sigma);


              mean     = fitResosEPhi[i][ebin][ebin2]->GetParameter(1);
              meanErr  = fitResosEPhi[i][ebin][ebin2]->GetParError(1);
              sigma    = fitResosEPhi[i][ebin][ebin2]->GetParameter(2);
              sigmaErr = fitResosEPhi[i][ebin][ebin2]->GetParError(2);
            }

            Double_t reso     = sigma/mean;
            Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
            histResoEPhiVsPhi[i][ebin]->SetBinContent(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[ebin2]),  100*reso);
            histResoEPhiVsPhi[i][ebin]->SetBinError(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[ebin2]), resoErr);

            histMeanEPhiVsPhi[i][ebin]->SetBinContent(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[ebin2]), mean);
            histMeanEPhiVsPhi[i][ebin]->SetBinError(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[ebin2]), meanErr);
            histSigmaEPhiVsPhi[i][ebin]->SetBinContent(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[ebin2]), sigma);
            histSigmaEPhiVsPhi[i][ebin]->SetBinError(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[ebin2]), sigmaErr);  

          }  
          
        }
      }
    }
  }

  
  for (Int_t i = 0; i < nParticles-1; i++){

    if (!useParEnetPhi[i]) continue;

     TCanvas* phidistE = new TCanvas("phidistE","",0,0,1100,1000);
     DrawGammaCanvasSettings( phidistE, 0.095, 0.01, 0.01, 0.105);
     phidistE->Divide(1,2);
     phidistE->cd(1);
     
     Int_t pC = 0; 
     TLegend* legendEPhi  = GetAndSetLegend2(0.2,0.2, 0.35, 0.9,0.4*textSizeLabelsPixel, 1, "", 43, 0.2);
     for (int ebin = 0; ebin < nEne; ebin++){
        if (!useEnePhi[i][ebin]) continue;
        
        drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsPixel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        
        if(useEnePhi[i][ebin]){

          SetStyleHistoTH1ForGraphs(histMeanEPhiVsPhi[i][ebin], "#phi ^{MC}","(#phi ^{rec} - #phi ^{MC})_{Mean}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
          histMeanEPhiVsPhi[i][ebin]->SetLineColor(pC+1);
          histMeanEPhiVsPhi[i][ebin]->GetYaxis()->SetDecimals(kTRUE);
          histMeanEPhiVsPhi[i][ebin]->GetYaxis()->SetMaxDigits(5);
          if(pC==0) {histMeanEPhiVsPhi[i][ebin]->Draw("");}
          else  {histMeanEPhiVsPhi[i][ebin]->Draw("SAME");}
          pC++;

          legendEPhi->AddEntry(histMeanEPhiVsPhi[i][ebin], Form("E^{MC} =  %1.1f", partE[ebin]),"lep");
        }   
     }
     legendEPhi->Draw();
   

     phidistE->cd(2);
     
     pC = 0; 
     TLegend* legendEEtaSigma  = GetAndSetLegend2(0.2,0.15, 0.35, 0.95,0.4*textSizeLabelsPixel, 1, "", 43, 0.2);
     for (int ebin = 0; ebin < nEne; ebin++){
        if (!useEnePhi[i][ebin]) continue;
        if(useEnePhi[i][ebin]){
          SetStyleHistoTH1ForGraphs(histSigmaEPhiVsPhi[i][ebin], "#phi ^{MC}","(#phi ^{rec} - #phi ^{MC})_{width}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
          histSigmaEPhiVsPhi[i][ebin]->SetLineColor(pC+1);
          histSigmaEPhiVsPhi[i][ebin]->GetYaxis()->SetDecimals(kTRUE);
          histSigmaEPhiVsPhi[i][ebin]->GetYaxis()->SetMaxDigits(5);
          if(pC==0) {histSigmaEPhiVsPhi[i][ebin]->Draw("");}
          else  {histSigmaEPhiVsPhi[i][ebin]->Draw("SAME");}
          pC++;

          legendEEtaSigma->AddEntry(histSigmaEPhiVsPhi[i][ebin], Form("E^{MC} =  %1.1f", partE[ebin]),"lep");
         }   
      }
      legendEEtaSigma->Draw();
      phidistE->Print(Form("%s/%sphivsE_%s.%s", outputDir.Data(), readNames[i].Data(),caloNameRead.Data(),suffix.Data()));
  }

  
  for (Int_t i = 0; i < nParticles; i++){
    Int_t nbins = 100; 
    Double_t mineta = -5.; Double_t maxeta = 5.; 
    histResoEtaVsEta[i]       = new TH1F(Form("h_%sreso_Eta", readNames[i].Data()), "", nbins, mineta, maxeta);
    histMeanEtaVsEta[i]       = new TH1F(Form("h_%sreso_EtaMean", readNames[i].Data()), "", nbins, mineta, maxeta);
    histSigmaEtaVsEta[i]       = new TH1F(Form("h_%sreso_EtaSigma", readNames[i].Data()), "", nbins, mineta, maxeta);

     Double_t minphi = -4.; Double_t maxphi = 4.; 
    // ======== Phi =================================//
    histResoPhiVsPhi[i]       = new TH1F(Form("h_%sreso_Phi", readNames[i].Data()), "", nbins, minphi, maxphi);
    histMeanPhiVsPhi[i]       = new TH1F(Form("h_%sreso_PhiMean", readNames[i].Data()), "", nbins, minphi, maxphi);
    histSigmaPhiVsPhi[i]       = new TH1F(Form("h_%sreso_PhiSigma", readNames[i].Data()), "", nbins, minphi, maxphi);

    TString tempName = Form("%s/h_ERH_Reso_%sEta_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempName = Form("%s/h_RH_ResoComb_%sEta_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );

    hist2DResoEta[i]    = (TH2F*)inputFile->Get(tempName.Data());

    // ======== Phi =================================//

    TString tempNamePhi = Form("%s/h_ERH_Reso_%sPhi_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempNamePhi = Form("%s/h_RH_ResoComb_%sPhi_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );

    hist2DResoPhi[i]    = (TH2F*)inputFile->Get(tempNamePhi.Data());
    if (hist2DResoPhi[i]->GetEntries() > 0){

     TCanvas* phidist = new TCanvas("phidist","",0,0,1100,1000);
     DrawGammaCanvasSettings( phidist, 0.095, 0.01, 0.01, 0.105);
     SetStyleHistoTH2ForGraphs(hist2DResoPhi[i], "#phi ^{MC}","(#phi ^{rec} - #phi ^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
     hist2DResoPhi[i]->GetYaxis()->SetNdivisions(505,kTRUE);
     hist2DResoPhi[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
     hist2DResoPhi[i]->Draw("colz");
     phidist->Print(Form("%s/%sphi_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

     for (int ebin = 0; ebin < nPhi; ebin++){

      hist1DResoPhibins[i][ebin]    = (TH1F*)hist2DResoPhi[i]->ProjectionY(Form("projection_dPhi_%s%d", readNames[i].Data(), ebin), 
                                                                        hist2DResoPhi[i]->GetXaxis()->FindBin(partPhi[ebin]-0.125), hist2DResoPhi[i]->GetXaxis()->FindBin(partPhi[ebin]+0.125),"e");

      if (hist1DResoPhibins[i][ebin]->GetEntries() > 10 ){

          usePartPhi[i]           = 1;
          usePhi[i][ebin]         = 1;
          nActivePhipart[i]++;
         
          Double_t mean     = FindLargestBin1DHist(hist1DResoPhibins[i][ebin], -2, 2);
          Double_t meanErr  = hist1DResoPhibins[i][ebin]->GetMeanError();
          Double_t sigma    = hist1DResoPhibins[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoPhibins[i][ebin]->GetRMSError();

          if (doFitting){
                
            fitResosPhi[i][ebin]     = new TF1(Form("fit_resophi_%s%d", readNames[i].Data(), ebin), "gaus", mean-3*sigma, mean+2*sigma);
            fitResosPhi[i][ebin]->SetParLimits(1, 0.6*mean, 1.2*mean);
            fitResosPhi[i][ebin]->SetParameter(1, mean);
            fitResosPhi[i][ebin]->SetParLimits(2, 0.6*sigma,sigma);
            fitResosPhi[i][ebin]->SetParameter(2, sigma);
            hist1DResoPhibins[i][ebin]->Fit(fitResosPhi[i][ebin],"L0RMEQ","",mean-1*sigma, mean+1.2*sigma);

            mean     = fitResosPhi[i][ebin]->GetParameter(1);
            meanErr  = fitResosPhi[i][ebin]->GetParError(1);
            sigma    = fitResosPhi[i][ebin]->GetParameter(2);
            sigmaErr = fitResosPhi[i][ebin]->GetParError(2);

          }
          
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          histResoPhiVsPhi[i]->SetBinContent(histResoPhiVsPhi[i]->FindBin(partPhi[ebin]),  100*reso);
          histResoPhiVsPhi[i]->SetBinError(histResoPhiVsPhi[i]->FindBin(partPhi[ebin]), resoErr);

          histMeanPhiVsPhi[i]->SetBinContent(histResoPhiVsPhi[i]->FindBin(partPhi[ebin]), mean);
          histMeanPhiVsPhi[i]->SetBinError(histResoPhiVsPhi[i]->FindBin(partPhi[ebin]), meanErr);
          histSigmaPhiVsPhi[i]->SetBinContent(histResoPhiVsPhi[i]->FindBin(partPhi[ebin]), sigma);
          histSigmaPhiVsPhi[i]->SetBinError(histResoPhiVsPhi[i]->FindBin(partPhi[ebin]), sigmaErr);  
      }  
    }

     TCanvas* phivar = new TCanvas("phivar","",0,0,1100,1000);
     phivar->Divide(1,2);
     DrawGammaCanvasSettings( phivar, 0.095, 0.01, 0.01, 0.105);
     phivar->cd(1);
     SetStyleHistoTH1ForGraphs(histMeanPhiVsPhi[i], "#phi","(#phi ^{rec} - #phi ^{MC})_{mean}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.75);
     DrawGammaSetMarker(histMeanPhiVsPhi[i], 20, 1.5, kBlack, kBlack);    
     histMeanPhiVsPhi[i]->GetYaxis()->SetDecimals(kTRUE);
     histMeanPhiVsPhi[i]->GetYaxis()->SetMaxDigits(5);
     histMeanPhiVsPhi[i]->Draw("PE");
     phivar->cd(2);
     SetStyleHistoTH1ForGraphs(histSigmaPhiVsPhi[i], "#phi","(#phi ^{rec} - #phi ^{MC})_{width}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.75);
     DrawGammaSetMarker(histSigmaPhiVsPhi[i], 20, 1.5, kBlack, kBlack);    
     histSigmaPhiVsPhi[i]->GetYaxis()->SetDecimals(kTRUE);
     histSigmaPhiVsPhi[i]->GetYaxis()->SetMaxDigits(5);
     histSigmaPhiVsPhi[i]->Draw("PE");
     phivar->Print(Form("%s/%sphi-mean-width_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 


    }


    // ======== Eta =================================//

    if (hist2DResoEta[i]->GetEntries() > 0){

     TCanvas* etadist = new TCanvas("etadist","",0,0,1100,1000);
     DrawGammaCanvasSettings( etadist, 0.095, 0.01, 0.01, 0.105);
     SetStyleHistoTH2ForGraphs(hist2DResoEta[i], "#eta ^{MC}","(#eta ^{rec} - #eta ^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
     hist2DResoEta[i]->GetYaxis()->SetNdivisions(505,kTRUE);
     hist2DResoEta[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
     hist2DResoEta[i]->Draw("colz");
     etadist->Print(Form("%s/%seta_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

     for (int ebin = 0; ebin < nEta; ebin++){

      hist1DResoEtabins[i][ebin]    = (TH1F*)hist2DResoEta[i]->ProjectionY(Form("projection_dEta_%s%d", readNames[i].Data(), ebin), 
                                                                        hist2DResoEta[i]->GetXaxis()->FindBin(partEta[ebin]-0.125), hist2DResoEta[i]->GetXaxis()->FindBin(partEta[ebin]+0.125),"e");

       if (hist1DResoEtabins[i][ebin]->GetEntries() > 10 ){
          usePartEta[i]            = 1;
          useEta[i][ebin]         = 1;
          nActiveEtapart[i]++;
          
          Double_t mean     = FindLargestBin1DHist(hist1DResoEtabins[i][ebin], -2, 2);
          Double_t meanErr  = hist1DResoEtabins[i][ebin]->GetMeanError();
          Double_t sigma    = hist1DResoEtabins[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoEtabins[i][ebin]->GetRMSError();


          if (doFitting){
                
            fitResosEta[i][ebin]     = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "gaus", mean-3*sigma, mean+2*sigma);
            fitResosEta[i][ebin]->SetParLimits(1, 0.6*mean, 1.2*mean);
            fitResosEta[i][ebin]->SetParameter(1, mean);
            fitResosEta[i][ebin]->SetParLimits(2, 0.6*sigma,sigma);
            fitResosEta[i][ebin]->SetParameter(2, sigma);
            hist1DResoEtabins[i][ebin]->Fit(fitResosEta[i][ebin],"L0RMEQ","",mean-1*sigma, mean+1.2*sigma);

            mean     = fitResosEta[i][ebin]->GetParameter(1);
            meanErr  = fitResosEta[i][ebin]->GetParError(1);
            sigma    = fitResosEta[i][ebin]->GetParameter(2);
            sigmaErr = fitResosEta[i][ebin]->GetParError(2);

          }
          
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          histResoEtaVsEta[i]->SetBinContent(histResoEtaVsEta[i]->FindBin(partEta[ebin]),  100*reso);
          histResoEtaVsEta[i]->SetBinError(histResoEtaVsEta[i]->FindBin(partEta[ebin]), resoErr);
          histMeanEtaVsEta[i]->SetBinContent(histResoEtaVsEta[i]->FindBin(partEta[ebin]), mean);
          histMeanEtaVsEta[i]->SetBinError(histResoEtaVsEta[i]->FindBin(partEta[ebin]), meanErr);
          histSigmaEtaVsEta[i]->SetBinContent(histResoEtaVsEta[i]->FindBin(partEta[ebin]), sigma);
          histSigmaEtaVsEta[i]->SetBinError(histResoEtaVsEta[i]->FindBin(partEta[ebin]), sigmaErr);  

       }
     }
  }
}


  // ======== Eta =================================//


   for (Int_t i = 0; i < nParticles-1; i++){

    if (!usePartEta[i]) continue;

    TCanvas* cPNG;
    split_canvas(cPNG, Form("cPNG%d",i), nActiveEtapart[i]);

    Int_t padnum =0;

    for(Int_t ebin=0; ebin< nEta; ebin++){

      if (!useEta[i][ebin]) continue;
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
    
      SetStyleHistoTH1ForGraphs(hist1DResoEtabins[i][ebin], "(#eta ^{rec} - #eta ^{MC})", "counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);

      hist1DResoEtabins[i][ebin]->GetXaxis()->SetRangeUser(-0.3,0.3);
      hist1DResoEtabins[i][ebin]->GetYaxis()->SetRangeUser(0.8,hist1DResoEtabins[i][ebin]->GetMaximum()*2);
      
      DrawGammaSetMarker(hist1DResoEtabins[i][ebin], 20, 1.5, kBlack, kBlack);      
      hist1DResoEtabins[i][ebin]->Draw("p,e");

      if (fitResosEta[i][ebin] && doFitting){
        DrawGammaSetMarkerTF1( fitResosEta[i][ebin], 7, 2, kBlue+1);
        fitResosEta[i][ebin]->Draw("same");
      }
    
      
      if (padnum == 0){
        drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.1f< #eta< %1.1f",minEtaMax,maxEtaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
       
      drawLatexAdd(Form("%0.3f < #eta < %0.3f ",partEta[ebin]-0.125, partEta[ebin]+0.125),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
     
      Int_t bin = histResoEtaVsEta[i]->FindBin(partEta[ebin]);
      DrawGammaLines(histResoEtaVsEta[i]->GetBinContent(bin), 
                     histResoEtaVsEta[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoEtabins[i][ebin]->GetMaximum(), 1, kGray+2, 2);
      DrawGammaLines(histMeanEtaVsEta[i]->GetBinContent(bin)-histMeanEtaVsEta[i]->GetBinContent(bin), 
                     histMeanEtaVsEta[i]->GetBinContent(bin)-histMeanEtaVsEta[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoEtabins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(histSigmaEtaVsEta[i]->GetBinContent(bin)+histSigmaEtaVsEta[i]->GetBinContent(bin), 
                     histSigmaEtaVsEta[i]->GetBinContent(bin)+histSigmaEtaVsEta[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoEtabins[i][ebin]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNG->Print(Form("%s/%setadist_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));
  }

  // ======== Phi =================================//

  for (Int_t i = 0; i < nParticles-1; i++){
    if (!usePartPhi[i]) continue;
    TCanvas* cPNGPhi;
    split_canvas(cPNGPhi, Form("cPNGPhi%d",i), nActivePhipart[i]);
    Int_t padnum =0;
      
    for(Int_t ebin=0; ebin< nPhi; ebin++){

      if (!usePhi[i][ebin]) continue;
      cPNGPhi->cd(padnum+1);
      DrawVirtualPadSettings( cPNGPhi->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
    
      SetStyleHistoTH1ForGraphs(hist1DResoPhibins[i][ebin], "(#phi ^{rec} - #phi ^{MC})", "counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);

      hist1DResoPhibins[i][ebin]->GetXaxis()->SetRangeUser(-0.3,0.3);
      hist1DResoPhibins[i][ebin]->GetYaxis()->SetRangeUser(0.8,hist1DResoPhibins[i][ebin]->GetMaximum()*2);

      DrawGammaSetMarker(hist1DResoPhibins[i][ebin], 20, 1.5, kBlack, kBlack);      
      hist1DResoPhibins[i][ebin]->Draw("p,e");

      if (fitResosPhi[i][ebin] && doFitting){
        DrawGammaSetMarkerTF1( fitResosPhi[i][ebin], 7, 2, kBlue+1);
        fitResosPhi[i][ebin]->Draw("same");
      }
    
      
      if (padnum == 0){
        drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.1f< #phi< %1.1f",minPhiMax,maxPhiMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
       
      drawLatexAdd(Form("%0.3f < #phi < %0.3f ",partPhi[ebin]-0.125, partPhi[ebin]+0.125),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
     
      Int_t bin = histResoPhiVsPhi[i]->FindBin(partPhi[ebin]);
      DrawGammaLines(histResoPhiVsPhi[i]->GetBinContent(bin), 
                     histResoEtaVsEta[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoPhibins[i][ebin]->GetMaximum(), 1, kGray+2, 2);
      DrawGammaLines(histMeanPhiVsPhi[i]->GetBinContent(bin)-histMeanPhiVsPhi[i]->GetBinContent(bin), 
                     histMeanPhiVsPhi[i]->GetBinContent(bin)-histMeanPhiVsPhi[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoPhibins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(histSigmaPhiVsPhi[i]->GetBinContent(bin)+histSigmaPhiVsPhi[i]->GetBinContent(bin), 
                     histSigmaPhiVsPhi[i]->GetBinContent(bin)+histSigmaPhiVsPhi[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoPhibins[i][ebin]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNGPhi->Print(Form("%s/%sphidist_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));
  }



}
