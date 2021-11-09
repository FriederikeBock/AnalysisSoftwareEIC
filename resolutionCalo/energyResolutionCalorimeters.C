#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

Double_t linTermECal[3][2]    = {{1,3}, {1,3}, {1,3}};
Double_t linTermHCal[3][2]    = {{6,10}, {7,10}, {7,10}};
Double_t sqrtETermECal[3][2]  = {{2,3}, {7,10}, {7,10}};
Double_t sqrtETermHCal[3][2]  = {{45,50}, {85,100}, {35,50}};

Double_t maxResPerEM[3]       = {12.5, 17.5, 20.5};
Double_t maxResPerHad[3]      = {62, 82, 72};
Double_t maxResEM[3]          = {1.125, 1.175, 1.205};
Double_t maxResHad[3]         = {1.25, 1.85, 1.65};


Double_t resFit1oEEM[3][2]    = {{1/TMath::Sqrt(20.), 1.45}, {-1, -1}, {-1, -1} };
Double_t resFitEM[3][2]       = {{0.2, 20}, {-1, -1}, {-1, -1} };
Double_t meanFitEM[3][2]      = {{10,-1}, { 10, -1}, {10, 20} };
Double_t meanFitHad[3][2]     = {{10,-1}, { 10, -1}, {10, 20} };

Double_t minEnergy[3][2]      = {{0,2.2}, { 0, 2.5}, {0, 1.5} };


void energyResolutionCalorimeters(
    TString inputFileName     = "",
    TString caloName          = "",
    TString clusterizerName   = "",
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
  Double_t maxReso    = 1.25;
  TString detLabel    = "";
  Bool_t isEMCal      = 0;
  Int_t rebinRes      = 1;
  TString caloNameRead = caloName;
  Bool_t isComb         = 0;
  Int_t dirCal          = 0;  // 0- bwd, 1- cent, 2-fwd
  Int_t exEbin          = 15;
  Int_t exEtabin        = 2;
    
  if (caloName.BeginsWith("EEMC") ){
    detLabel      = "EEMC (PbWO_{4} crystal)";
    caloNameRead  = "EEMC";
    minEtaMax     = -4;
    maxEtaMax     = -1.7;
    isEMCal       = 1;
    exEtabin      = 2;
    if (caloName.CompareTo("EEMC-wMat") == 0){
      detLabel      = "EEMC, w/ mat. infront";
    }
  } else if (caloName.BeginsWith("FEMC") ){
    detLabel  = "FEMC (PHENIX re-use)";
    caloNameRead = "FEMC";
    minEtaMax     = 1.4;
    maxEtaMax     = 4;
    isEMCal       = 1;
    dirCal        = 2;
    exEtabin      = 11;
    if (caloName.CompareTo("FEMC-wMat") == 0){
      detLabel  = "FEMC, w/ mat. infront";
    }
  } else if (caloName.BeginsWith("BECAL") ){
    detLabel      = "BECAL (Sci-glass)";
    caloNameRead  = "BECAL";
    minEtaMax     = -1.8;
    maxEtaMax     = 1.3;
    isEMCal       = 1;
    exEtabin      = 7;
    dirCal        = 1;
    if (caloName.CompareTo("BECAL-wMat") == 0){
      detLabel  = "BECAL, w/ mat. infront";
    }
  } else if (caloName.BeginsWith("LFHCAL") ){
    detLabel      = "LFHCAL";
    caloNameRead  = "LFHCAL";
    minEtaMax     = 1.1;
    maxEtaMax     = 4;
    rebinRes      = 2;
    exEtabin      = 11;
    dirCal        = 2;
    if (caloName.CompareTo("LFHCAL-wMat") == 0){
      detLabel    = "LFHCAL, FEMC infront";
    } else if (caloName.CompareTo("LFHCAL-FEMC") == 0){
      detLabel    = "LFHCAL+FEMC";
      isComb      = 1;
    }
  } else if (caloName.BeginsWith("EHCAL") ){
    detLabel      = "EHCAL (STAR re-use)";
    caloNameRead  = "EHCAL";
    minEtaMax     = -4;
    maxEtaMax     = -1.7;
    exEtabin      = 2;
    if (caloName.CompareTo("EHCAL-wMat") == 0){
      detLabel  = "EHCAL, EEMC infront";
    } else if (caloName.CompareTo("EHCAL-EEMC") == 0){
      detLabel  = "EHCAL+EEMC";
      isComb      = 1;
    }
  } else if (caloName.BeginsWith("CHCAL") ){
    detLabel      = "oHCAL, iHCal infront";
    caloNameRead  = "HCALOUT";
    minEtaMax     = -1.8;
    maxEtaMax     = 1.2;
    exEtabin      = 2;
    dirCal        = 1;
    if (caloName.CompareTo("CHCAL-comb") == 0){
      detLabel    = "oHCAL+iHCal+BECAL";
      isComb      = 1;
    }
  } else if (caloName.BeginsWith("HCALIN") ){
    detLabel      = "iHCal";
    caloNameRead  = "HCALIN";
    minEtaMax     = -1.8;
    maxEtaMax     = 1.2;
    exEtabin      = 2;
    dirCal        = 1;
    if (caloName.CompareTo("HCALIN-wMat") == 0){
      detLabel  = "iHCal-w/ mat infront";
    }
  } else {
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
  }
  // determine eta bins
  Int_t binEtaMin = 0;
  Int_t binEtaMax = 15;
  
  TString outputDir                 = Form("plotsResolution%s",caloName.Data());
  if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);


  //************************** Read data **************************************************
  const int nEne = 34;
  
  
  const static Double_t partE[]   = {   0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 
                                        5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 
                                        13, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 50, 
                                        67.5, 75, 100, 150
                                    };

  TGraphErrors* graphReqEM          = new TGraphErrors(nEne); 
  TGraphErrors* graphReq1oEEM       = new TGraphErrors(nEne); 
  TGraphErrors* graphReqHad         = new TGraphErrors(nEne); 
  TGraphErrors* graphReq1oEHad      = new TGraphErrors(nEne); 
  for (Int_t i= 0; i < nEne; i++){
//     cout << "ebin:" << partE[i] << endl;
    Double_t value[4]             = {0};
    value[0]  = sqrtETermECal[dirCal][0]/TMath::Sqrt(partE[i])+linTermECal[dirCal][0];
    value[1]  = sqrtETermECal[dirCal][1]/TMath::Sqrt(partE[i])+linTermECal[dirCal][0];
    value[2]  = sqrtETermECal[dirCal][0]/TMath::Sqrt(partE[i])+linTermECal[dirCal][1];
    value[3]  = sqrtETermECal[dirCal][1]/TMath::Sqrt(partE[i])+linTermECal[dirCal][1];
    Double_t valueHad[4]             = {0};
    valueHad[0]  = sqrtETermHCal[dirCal][0]/TMath::Sqrt(partE[i])+linTermHCal[dirCal][0];
    valueHad[1]  = sqrtETermHCal[dirCal][1]/TMath::Sqrt(partE[i])+linTermHCal[dirCal][0];
    valueHad[2]  = sqrtETermHCal[dirCal][0]/TMath::Sqrt(partE[i])+linTermHCal[dirCal][1];
    valueHad[3]  = sqrtETermHCal[dirCal][1]/TMath::Sqrt(partE[i])+linTermHCal[dirCal][1];
//     cout << value[0] << "\t" << value[1] << "\t"<< value[2] << "\t"<< value[3] << "\t"<< endl;
    Double_t max = -10000;
    Double_t min = 10000;
    for (Int_t k = 0; k < 4; k++){ 
        if (min > value[k]) min = value[k];
        if (max < value[k]) max = value[k];
    }
    graphReqEM->SetPoint(i, partE[i], (min+max)/2 );
    graphReqEM->SetPointError(i, 0.05, TMath::Abs(min-max)/2 );
    graphReq1oEEM->SetPoint(i, 1/TMath::Sqrt(partE[i]), (min+max)/2 );
    graphReq1oEEM->SetPointError(i, 0.005, TMath::Abs(min-max)/2 );
    
    max = -10000;
    min = 10000;
    for (Int_t k = 0; k < 4; k++){ 
        if (min > valueHad[k]) min = valueHad[k];
        if (max < valueHad[k]) max = valueHad[k];
    }
    graphReqHad->SetPoint(i, partE[i], (min+max)/2 );
    graphReqHad->SetPointError(i, 0.05, TMath::Abs(min-max)/2 );
    graphReq1oEHad->SetPoint(i, 1/TMath::Sqrt(partE[i]), (min+max)/2 );
    graphReq1oEHad->SetPointError(i, 0.005, TMath::Abs(min-max)/2 );
  }
  graphReq1oEEM->Sort();
  DrawGammaSetMarkerTGraphErr(  graphReqEM, 0, 0, kBlue-10, kBlue-10, 1, kTRUE, kBlue-10, kFALSE);
  DrawGammaSetMarkerTGraphErr(  graphReq1oEEM, 0, 0, kBlue-10, kBlue-10, 1, kTRUE, kBlue-10, kFALSE);
  DrawGammaSetMarkerTGraphErr(  graphReqHad, 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
  DrawGammaSetMarkerTGraphErr(  graphReq1oEHad, 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
  
  Int_t useE[7][nEne]             = {{ 0}};
  Int_t usePart[7]                = { 0};
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
  Int_t nActiveEpart[7]             = {0};
  Int_t nActiveEtapart[7]           = {0};
  
  TFile* inputFile                  = new TFile(inputFileName.Data());
  TH2D* hist2DResoE[7]              = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH2D* hist2DResoEhighest[7]       = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1D* hist1DResoEbins[7][nEne]    = {{NULL}};
  TH1D* hist1DResoEbinshighests[7][nEne]    = {{NULL}};
  TH1F* histSigmaEVsE[7]            = {NULL};
  TH1F* histFWHMEVsE[7]             = {NULL};
  TH1F* histResoEVsE[7]             = {NULL};
  TH1F* histResoFEVsE[7]            = {NULL};
  TH1F* histMeanEVsE[7]             = {NULL};
  TH1F* histResoEVs1oE[7]           = {NULL};
  TH1F* histResoFEVs1oE[7]          = {NULL};
  TF1* fitMeanE[7]                  = {NULL};
  TF1* fitResoFE[7]                  = {NULL};
  TF1* fitResoF1oE[7]                = {NULL};
  TF1* fitResoE[7]                  = {NULL};
  TF1* fitReso1oE[7]                = {NULL};
  TF1* fitResoEbins[7][nEne]        = {{NULL}};
  TH1F* histSigmaEVsEhighest[7]            = {NULL};
  TH1F* histFWHMEVsEhighest[7]             = {NULL};
  TH1F* histResoEVsEhighest[7]             = {NULL};
  TH1F* histResoFEVsEhighest[7]            = {NULL};
  TH1F* histMeanEVsEhighest[7]             = {NULL};
  TH1F* histResoEVs1oEhighest[7]           = {NULL};
  TH1F* histResoFEVs1oEhighest[7]          = {NULL};
  TF1* fitResoEhighest[7]                  = {NULL};
  TF1* fitReso1oEhighest[7]                = {NULL};
  TF1* fitResoFEhighest[7]                 = {NULL};
  TF1* fitResoF1oEhighest[7]               = {NULL};
  TF1* fitResoEbinshighest[7][nEne]        = {{NULL}};

  TH2D* hist2DResoEEta[7][nEta]              = {{NULL}};
  TH1D* hist1DResoEbinsEta[7][nEne][nEta]    = {{{NULL}}};
  TH1F* histSigmaEVsEEta[7][nEta]            = {{NULL}};
  TH1F* histFWHMEVsEEta[7][nEta]             = {{NULL}};
  TH1F* histResoEVsEEta[7][nEta]             = {{NULL}};
  TH1F* histResoFEVsEEta[7][nEta]            = {{NULL}};
  TH1F* histMeanEVsEEta[7][nEta]             = {{NULL}};
  TH1F* histResoEVs1oEEta[7][nEta]           = {{NULL}};
  TH1F* histResoFEVs1oEEta[7][nEta]          = {{NULL}};
  TF1* fitResoEEta[7][nEta]                  = {{NULL}};
  TF1* fitReso1oEEta[7][nEta]                = {{NULL}};
  TF1* fitResoFEEta[7][nEta]                  = {{NULL}};
  TF1* fitResoF1oEEta[7][nEta]                = {{NULL}};
  TF1* fitResoEbinsEta[7][nEne][nEta]        = {{NULL}};
  Int_t useEEta[7][nEne][nEta]               = {{{ 0}}};
  Int_t useEta[7][nEta]                      = {{0}};
  
  for (Int_t i = 0; i < nParticles; i++){
    histSigmaEVsE[i]      = new TH1F(Form("h_%ssigma_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histFWHMEVsE[i]       = new TH1F(Form("h_%sfwhm_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histMeanEVsE[i]       = new TH1F(Form("h_%smean_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVsE[i]       = new TH1F(Form("h_%sreso_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVs1oE[i]     = new TH1F(Form("h_%sreso_1oE", readNames[i].Data()), "", 1000, 0., 4.);
    histResoFEVsE[i]       = new TH1F(Form("h_%sresof_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoFEVs1oE[i]     = new TH1F(Form("h_%sresof_1oE", readNames[i].Data()), "", 1000, 0., 4.);
    histSigmaEVsEhighest[i]      = new TH1F(Form("h_%ssigma_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histFWHMEVsEhighest[i]       = new TH1F(Form("h_%sfwhm_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histMeanEVsEhighest[i]       = new TH1F(Form("h_%smean_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVsEhighest[i]       = new TH1F(Form("h_%sreso_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVs1oEhighest[i]     = new TH1F(Form("h_%sreso_1oEhighest", readNames[i].Data()), "", 1000, 0., 4.);
    histResoFEVsEhighest[i]       = new TH1F(Form("h_%sresof_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoFEVs1oEhighest[i]     = new TH1F(Form("h_%sresof_1oEhighest", readNames[i].Data()), "", 1000, 0., 4.);
    
    
    TString tempName = Form("%s/h_CRH_EReso_%sE_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempName = Form("%s/h_CRH_EResoComb_%sE_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    cout << tempName << endl;
    TString tempName3 = Form("%s/h_CRH_EReso_%sEhighest_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempName3 = tempName;
    cout << tempName3 << endl;
    hist2DResoE[i]            = (TH2D*)inputFile->Get(tempName.Data());
    hist2DResoEhighest[i]     = (TH2D*)inputFile->Get(tempName3.Data());
    if (!hist2DResoE[i]) continue;
    if (hist2DResoE[i]->GetEntries() > 0){
      hist2DResoE[i]->Sumw2();
      hist2DResoEhighest[i]->Sumw2();
      for (int ebin = 0; ebin < nEne; ebin++){
        hist1DResoEbins[i][ebin]    = (TH1D*)hist2DResoE[i]->ProjectionY(Form("projection_dE_%s%d", readNames[i].Data(), ebin), 
                                                                        hist2DResoE[i]->GetXaxis()->FindBin(partE[ebin]-0.05), hist2DResoE[i]->GetXaxis()->FindBin(partE[ebin]+0.05),"e"); 
        if (i < 2 && partE[ebin] < minEnergy[dirCal][0]) continue; // switch of low E bins for ECals
        if (i > 1 && partE[ebin] < minEnergy[dirCal][1]) continue; // switch of low E bins for HCals
        
        if (hist1DResoEbins[i][ebin]->GetEntries() > 10 ){
          useE[i][ebin]         = 1;
          usePart[i]            = 1;
          nActiveEpart[i]++;
          
          Double_t mean     = FindLargestBin1DHist(hist1DResoEbins[i][ebin], 0.4, 1.9);
          Double_t meanErr  = hist1DResoEbins[i][ebin]->GetMeanError();
//           Double_t sigma    = hist1DResoEbins[i][ebin]->GetStdDev();
//           Double_t sigmaErr = hist1DResoEbins[i][ebin]->GetStdDevError();
          Double_t sigma    = hist1DResoEbins[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoEbins[i][ebin]->GetRMSError();
          Double_t fwhm    = hist1DResoEbins[i][ebin]->GetRMS();
          Double_t fwhmErr = hist1DResoEbins[i][ebin]->GetRMSError();
          
          if (doFitting){
            cout << labelPart[i] <<"\t" << partE[ebin] << "\t" << mean << "\t" << sigma<< endl;;
//             fitResoEbins[i][ebin]     = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))", mean-3*sigma, mean+2*sigma);
            hist1DResoEbins[i][ebin]->Rebin(rebinRes);
            if (i > 1){
              fitResoEbins[i][ebin]     = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "gaus", mean-3*sigma, mean+2*sigma);
              fitResoEbins[i][ebin]->SetParLimits(1, 0.6*mean, 1.2*mean);
              fitResoEbins[i][ebin]->SetParLimits(1, 0.6*mean, 1.2*mean);
              fitResoEbins[i][ebin]->SetParameter(1, mean);
              if (!isEMCal && dirCal == 1)
                fitResoEbins[i][ebin]->SetParLimits(2, 0.6*sigma,sigma);
              else
                fitResoEbins[i][ebin]->SetParLimits(2, 0.3*sigma,2*sigma);
              fitResoEbins[i][ebin]->SetParameter(2, sigma);
//             fitResoEbins[i][ebin]->SetParLimits(0, 0.9*hist1DResoEbins[i][ebin]->GetMaximum(), 10*hist1DResoEbins[i][ebin]->GetMaximum());
//             fitResoEbins[i][ebin]->SetParameter(0, hist1DResoEbins[i][ebin]->GetMaximum());
            } else  {
              fitResoEbins[i][ebin]    = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "crystalball", mean-3*sigma, mean+2*sigma);
              fitResoEbins[i][ebin]->SetParameters(1*hist1DResoEbins[i][ebin]->GetMaximum(),hist1DResoEbins[i][ebin]->GetMean(),0.5*hist1DResoEbins[i][ebin]->GetRMS(),0.2*hist1DResoEbins[i][ebin]->GetRMS(),1*hist1DResoEbins[i][ebin]->GetRMS());
             
            }
            hist1DResoEbins[i][ebin]->Fit(fitResoEbins[i][ebin],"L0RME","",mean-1*sigma, mean+1.2*sigma);
//             hist1DResoEbins[i][ebin]->Fit(fitResoEbins[i][ebin],"L0RME", "", fitResoEbins[i][ebin]->GetParameter(1)-0.5*fitResoEbins[i][ebin]->GetParameter(2), fitResoEbins[i][ebin]->GetParameter(1)+1*fitResoEbins[i][ebin]->GetParameter(2));
            
            mean     = fitResoEbins[i][ebin]->GetParameter(1);
            meanErr  = fitResoEbins[i][ebin]->GetParError(1);
            sigma    = fitResoEbins[i][ebin]->GetParameter(2);
            sigmaErr = fitResoEbins[i][ebin]->GetParError(2);
            if (i > 1){
              fwhm    = sigma;
              fwhmErr = sigmaErr;
            } else {
              fwhm  = (fitResoEbins[i][ebin]->GetParameter(1) - fitResoEbins[i][ebin]->GetX(0.5*fitResoEbins[i][ebin]->GetMaximum(), 0.2 , fitResoEbins[i][ebin]->GetParameter(1)));
              Double_t relEM = fitResoEbins[i][ebin]->GetParError(1)/ fitResoEbins[i][ebin]->GetParameter(1);
              Double_t relES = fitResoEbins[i][ebin]->GetParError(2)/ fitResoEbins[i][ebin]->GetParameter(2);
              Double_t relEL = fitResoEbins[i][ebin]->GetParError(3)/ fitResoEbins[i][ebin]->GetParameter(3);
              Double_t relEN = fitResoEbins[i][ebin]->GetParError(4)/ fitResoEbins[i][ebin]->GetParameter(4);
              fwhmErr = fwhm*TMath::Sqrt(relEM*relEM+relES*relES+relEL*relEL+relEN*relEN);
            }
          }
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          Double_t resoF    = fwhm/mean;
          Double_t resoFErr = TMath::Sqrt(TMath::Power(fwhmErr/fwhm,2)+TMath::Power(meanErr/mean,2));
          
          histResoEVsE[i]->SetBinContent(histResoEVsE[i]->FindBin(partE[ebin]),  100*reso);
          histResoEVsE[i]->SetBinError(histResoEVsE[i]->FindBin(partE[ebin]), resoErr);
          histResoEVs1oE[i]->SetBinContent(histResoEVs1oE[i]->FindBin(1./TMath::Sqrt(partE[ebin])),100*reso);
          histResoEVs1oE[i]->SetBinError(histResoEVs1oE[i]->FindBin(1./TMath::Sqrt(partE[ebin])),resoErr);
          histResoFEVsE[i]->SetBinContent(histResoFEVsE[i]->FindBin(partE[ebin]),  100*resoF);
          histResoFEVsE[i]->SetBinError(histResoFEVsE[i]->FindBin(partE[ebin]), resoFErr);
          histResoFEVs1oE[i]->SetBinContent(histResoFEVs1oE[i]->FindBin(1./TMath::Sqrt(partE[ebin])),100*resoF);
          histResoFEVs1oE[i]->SetBinError(histResoFEVs1oE[i]->FindBin(1./TMath::Sqrt(partE[ebin])),resoFErr);

          histMeanEVsE[i]->SetBinContent(histResoEVsE[i]->FindBin(partE[ebin]), mean);
          histMeanEVsE[i]->SetBinError(histResoEVsE[i]->FindBin(partE[ebin]), meanErr);
          histSigmaEVsE[i]->SetBinContent(histResoEVsE[i]->FindBin(partE[ebin]), sigma);
          histSigmaEVsE[i]->SetBinError(histResoEVsE[i]->FindBin(partE[ebin]), sigmaErr);
          histFWHMEVsE[i]->SetBinContent(histResoEVsE[i]->FindBin(partE[ebin]), fwhm);
          histFWHMEVsE[i]->SetBinError(histResoEVsE[i]->FindBin(partE[ebin]), fwhmErr);
        } else {
          useE[i][ebin] = 0;
        }
      }
      for (int ebin = 0; ebin < nEne; ebin++){
        hist1DResoEbinshighests[i][ebin]    = (TH1D*)hist2DResoEhighest[i]->ProjectionY(Form("projection_dEhighest_%s%d", readNames[i].Data(), ebin), 
                                                                        hist2DResoEhighest[i]->GetXaxis()->FindBin(partE[ebin]-0.05), hist2DResoEhighest[i]->GetXaxis()->FindBin(partE[ebin]+0.05),"e"); 
        if (i < 2 && partE[ebin] < minEnergy[dirCal][0]) continue; // switch of low E bins for ECals
        if (i > 1 && partE[ebin] < minEnergy[dirCal][1]) continue; // switch of low E bins for HCals

        if (hist1DResoEbinshighests[i][ebin]->GetEntries() > 10 ){
          Double_t mean     = FindLargestBin1DHist(hist1DResoEbinshighests[i][ebin], 0.4, 1.9);
//           Double_t mean     = hist1DResoEbins[i][ebin]->GetMean();
          Double_t meanErr  = hist1DResoEbinshighests[i][ebin]->GetMeanError();
//           Double_t sigma    = hist1DResoEbins[i][ebin]->GetStdDev();
//           Double_t sigmaErr = hist1DResoEbins[i][ebin]->GetStdDevError();
          Double_t sigma    = hist1DResoEbinshighests[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoEbinshighests[i][ebin]->GetRMSError();
          Double_t fwhm     = hist1DResoEbinshighests[i][ebin]->GetRMS();
          Double_t fwhmErr  = hist1DResoEbinshighests[i][ebin]->GetRMSError();

          if (doFitting){
            hist1DResoEbinshighests[i][ebin]->Rebin(rebinRes);
            if (i > 1){
//             cout << labelPart[i] <<"\t" << partE[ebin] << "\t" << mean << "\t" << sigma<< endl;;
//             fitResoEbins[i][ebin]     = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))", mean-3*sigma, mean+2*sigma);
              fitResoEbinshighest[i][ebin]     = new TF1(Form("fit_resohighest_%s%d", readNames[i].Data(), ebin), "gaus", mean-3*sigma, mean+2*sigma);
              if (!isEMCal && dirCal == 1)
                fitResoEbinshighest[i][ebin]->SetParLimits(1, 0.6*mean, 1.2*mean);
              else 
                fitResoEbinshighest[i][ebin]->SetParLimits(1, 0.4, 1.2);
              fitResoEbinshighest[i][ebin]->SetParameter(1, mean);
              fitResoEbinshighest[i][ebin]->SetParLimits(2, 0.6*sigma,sigma);
              fitResoEbinshighest[i][ebin]->SetParameter(2, sigma);
//             fitResoEbins[i][ebin]->SetParLimits(0, 0.9*hist1DResoEbins[i][ebin]->GetMaximum(), 10*hist1DResoEbins[i][ebin]->GetMaximum());
//             fitResoEbins[i][ebin]->SetParameter(0, hist1DResoEbins[i][ebin]->GetMaximum());
            } else  {
              fitResoEbinshighest[i][ebin]    = new TF1(Form("fit_resohighest_%s%d", readNames[i].Data(), ebin), "crystalball", mean-3*sigma, mean+2*sigma);
              fitResoEbinshighest[i][ebin]->SetParameters(1*hist1DResoEbinshighests[i][ebin]->GetMaximum(),hist1DResoEbinshighests[i][ebin]->GetMean(),0.5*hist1DResoEbinshighests[i][ebin]->GetRMS(),0.2*hist1DResoEbinshighests[i][ebin]->GetRMS(),1*hist1DResoEbinshighests[i][ebin]->GetRMS());

            }
            hist1DResoEbinshighests[i][ebin]->Fit(fitResoEbinshighest[i][ebin],"L0RMEQ","",mean-1*sigma, mean+1.2*sigma);
//             hist1DResoEbins[i][ebin]->Fit(fitResoEbins[i][ebin],"L0RME", "", fitResoEbins[i][ebin]->GetParameter(1)-0.5*fitResoEbins[i][ebin]->GetParameter(2), fitResoEbins[i][ebin]->GetParameter(1)+1*fitResoEbins[i][ebin]->GetParameter(2));

            mean     = fitResoEbinshighest[i][ebin]->GetParameter(1);
            meanErr  = fitResoEbinshighest[i][ebin]->GetParError(1);
            sigma    = fitResoEbinshighest[i][ebin]->GetParameter(2);
            sigmaErr = fitResoEbinshighest[i][ebin]->GetParError(2);
            if (i > 1){
              fwhm    = sigma;
              fwhmErr = sigmaErr;
            } else {
              fwhm  = (fitResoEbinshighest[i][ebin]->GetParameter(1) - fitResoEbinshighest[i][ebin]->GetX(0.5*fitResoEbinshighest[i][ebin]->GetMaximum(), 0.2 , fitResoEbinshighest[i][ebin]->GetParameter(1)));
              Double_t relEM = fitResoEbinshighest[i][ebin]->GetParError(1)/ fitResoEbinshighest[i][ebin]->GetParameter(1);
              Double_t relES = fitResoEbinshighest[i][ebin]->GetParError(2)/ fitResoEbinshighest[i][ebin]->GetParameter(2);
              Double_t relEL = fitResoEbinshighest[i][ebin]->GetParError(3)/ fitResoEbinshighest[i][ebin]->GetParameter(3);
              Double_t relEN = fitResoEbinshighest[i][ebin]->GetParError(4)/ fitResoEbinshighest[i][ebin]->GetParameter(4);
              fwhmErr = fwhm*TMath::Sqrt(relEM*relEM+relES*relES+relEL*relEL+relEN*relEN);              
            }
          }
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          Double_t resoF    = fwhm/mean;
          Double_t resoFErr = TMath::Sqrt(TMath::Power(fwhmErr/fwhm,2)+TMath::Power(meanErr/mean,2));

          
          histResoEVsEhighest[i]->SetBinContent(histResoEVsEhighest[i]->FindBin(partE[ebin]),  100*reso);
          histResoEVsEhighest[i]->SetBinError(histResoEVsEhighest[i]->FindBin(partE[ebin]), resoErr);
          histResoEVs1oEhighest[i]->SetBinContent(histResoEVs1oEhighest[i]->FindBin(1./TMath::Sqrt(partE[ebin])),100*reso);
          histResoEVs1oEhighest[i]->SetBinError(histResoEVs1oEhighest[i]->FindBin(1./TMath::Sqrt(partE[ebin])),resoErr);
          histResoFEVsEhighest[i]->SetBinContent(histResoFEVsEhighest[i]->FindBin(partE[ebin]),  100*resoF);
          histResoFEVsEhighest[i]->SetBinError(histResoFEVsEhighest[i]->FindBin(partE[ebin]), resoFErr);
          histResoFEVs1oEhighest[i]->SetBinContent(histResoFEVs1oEhighest[i]->FindBin(1./TMath::Sqrt(partE[ebin])),100*resoF);
          histResoFEVs1oEhighest[i]->SetBinError(histResoFEVs1oEhighest[i]->FindBin(1./TMath::Sqrt(partE[ebin])),resoFErr);

          histMeanEVsEhighest[i]->SetBinContent(histResoEVsEhighest[i]->FindBin(partE[ebin]), mean);
          histMeanEVsEhighest[i]->SetBinError(histResoEVsEhighest[i]->FindBin(partE[ebin]), meanErr);
          histSigmaEVsEhighest[i]->SetBinContent(histResoEVsEhighest[i]->FindBin(partE[ebin]), sigma);
          histSigmaEVsEhighest[i]->SetBinError(histResoEVsEhighest[i]->FindBin(partE[ebin]), sigmaErr);
          histFWHMEVsEhighest[i]->SetBinContent(histResoEVsEhighest[i]->FindBin(partE[ebin]), fwhm);
          histFWHMEVsEhighest[i]->SetBinError(histResoEVsEhighest[i]->FindBin(partE[ebin]), fwhmErr);
        } else {
          useE[i][ebin] = 0;
        }
      }
      for (Int_t etbin = 0; etbin < nEta; etbin++){
        histSigmaEVsEEta[i][etbin]      = new TH1F(Form("h_%ssigma_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histFWHMEVsEEta[i][etbin]      = new TH1F(Form("h_%sfwhm_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histMeanEVsEEta[i][etbin]       = new TH1F(Form("h_%smean_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histResoEVsEEta[i][etbin]       = new TH1F(Form("h_%sreso_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histResoEVs1oEEta[i][etbin]     = new TH1F(Form("h_%sreso_1oE_Eta_%d", readNames[i].Data(), etbin), "", 1000, 0., 4.);
        histResoFEVsEEta[i][etbin]       = new TH1F(Form("h_%sresof_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histResoFEVs1oEEta[i][etbin]     = new TH1F(Form("h_%sresof_1oE_Eta_%d", readNames[i].Data(), etbin), "", 1000, 0., 4.);
        // loading higest energy hist for particle
        TString tempName2 = Form("%s/h_CRH_EhighestReso_%sE_Eta_%s_%s_%d",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data(), etbin );
        if (isComb) tempName2 = Form("%s/h_CRH_EResoComb_%sE_Eta_%s_%s_%d",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data(), etbin );
//         cout << tempName2 << endl;

        hist2DResoEEta[i][etbin]    = (TH2D*)inputFile->Get(tempName2.Data());
        if (!hist2DResoEEta[i][etbin]) continue;

        if (hist2DResoEEta[i][etbin]->GetEntries() > 0){
          hist2DResoEEta[i][etbin]->Sumw2();
          for (int ebin = 0; ebin < nEne; ebin++){
            
            hist1DResoEbinsEta[i][ebin][etbin]    = (TH1D*)hist2DResoEEta[i][etbin]->ProjectionY(Form("projection_dE_%s%d_%d", readNames[i].Data(), ebin, etbin), 
                                                                            hist2DResoEEta[i][etbin]->GetXaxis()->FindBin(partE[ebin]-0.05), hist2DResoEEta[i][etbin]->GetXaxis()->FindBin(partE[ebin]+0.05),"e"); 
            if (i < 2 && partE[ebin] < minEnergy[dirCal][0]) continue; // switch of low E bins for ECals
            if (i > 1 && partE[ebin] < minEnergy[dirCal][1]) continue; // switch of low E bins for HCals
            if (hist1DResoEbinsEta[i][ebin][etbin]->GetEntries() > 10){
              useEEta[i][ebin][etbin]  = 1;
              if (!useEta[i][etbin]) nActiveEtapart[i]++;
              useEta[i][etbin]         = 1;
              
              Double_t meanEta     = FindLargestBin1DHist(hist1DResoEbinsEta[i][ebin][etbin], 0.5, 1.9);
//               Double_t meanEta     = hist1DResoEbinsEta[i][ebin][etbin]->GetMean();
              Double_t meanEtaErr  = hist1DResoEbinsEta[i][ebin][etbin]->GetMeanError();
//               Double_t sigmaEta    = hist1DResoEbinsEta[i][ebin][etbin]->GetStdDev();
//               Double_t sigmaEtaErr = hist1DResoEbinsEta[i][ebin][etbin]->GetStdDevError();
              Double_t sigmaEta     = hist1DResoEbinsEta[i][ebin][etbin]->GetRMS();
              Double_t sigmaEtaErr  = hist1DResoEbinsEta[i][ebin][etbin]->GetRMSError();
              Double_t fwhmEta      = hist1DResoEbinsEta[i][ebin][etbin]->GetRMS();
              Double_t fwhmEtaErr   = hist1DResoEbinsEta[i][ebin][etbin]->GetRMSError();

              if (doFitting){
                hist1DResoEbinsEta[i][ebin][etbin]->Rebin(rebinRes);
                if (i > 1){
                  fitResoEbinsEta[i][ebin][etbin]     = new TF1(Form("fit_reso_%s%d_%d", readNames[i].Data(), ebin, etbin), "gaus", meanEta-3*sigmaEta, meanEta+2*sigmaEta);
  //                 fitResoEbinsEta[i][ebin][etbin]     = new TF1(Form("fit_reso_%s%d_%d", readNames[i].Data(), ebin, etbin), "(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))", meanEta-3*sigmaEta, meanEta+2*sigmaEta);
                  
                  fitResoEbinsEta[i][ebin][etbin]->SetParLimits(1, 0.6*meanEta, 1.2*meanEta);
                  fitResoEbinsEta[i][ebin][etbin]->SetParameter(1, meanEta);
                  fitResoEbinsEta[i][ebin][etbin]->SetParLimits(2, 0.7*sigmaEta, sigmaEta);
                  fitResoEbinsEta[i][ebin][etbin]->SetParameter(2, sigmaEta);
  //                 fitResoEbinsEta[i][ebin][etbin]->SetParLimits(0, 0.9*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 10*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum());
  //                 fitResoEbinsEta[i][ebin][etbin]->SetParameter(0, hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum());
                  fitResoEbinsEta[i][ebin][etbin]->SetParLimits(1, 0.6*fitResoEbins[i][ebin]->GetParameter(1), 1.4*fitResoEbins[i][ebin]->GetParameter(1));
                  fitResoEbinsEta[i][ebin][etbin]->SetParameter(1, fitResoEbins[i][ebin]->GetParameter(1));
                  fitResoEbinsEta[i][ebin][etbin]->SetParLimits(2, 0.6*fitResoEbins[i][ebin]->GetParameter(2), 1.4*fitResoEbins[i][ebin]->GetParameter(2));
                  fitResoEbinsEta[i][ebin][etbin]->SetParameter(2, fitResoEbins[i][ebin]->GetParameter(2));                
                } else  {
                  fitResoEbinsEta[i][ebin][etbin]    = new TF1(Form("fit_reso_%s%d_%d", readNames[i].Data(), ebin, etbin), "crystalball", meanEta-3*sigmaEta, meanEta+2*sigmaEta);
                  fitResoEbinsEta[i][ebin][etbin]->SetParameters(1*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(),hist1DResoEbinsEta[i][ebin][etbin]->GetMean(),0.5*hist1DResoEbinsEta[i][ebin][etbin]->GetRMS(),0.2*hist1DResoEbinsEta[i][ebin][etbin]->GetRMS(),1*hist1DResoEbinshighests[i][ebin]->GetRMS());
                  
                }
                hist1DResoEbinsEta[i][ebin][etbin]->Fit(fitResoEbinsEta[i][ebin][etbin],"L0RMEQ","", meanEta-1*sigmaEta, meanEta+1.2*sigmaEta);
//                 hist1DResoEbinsEta[i][ebin][etbin]->Fit(fitResoEbinsEta[i][ebin][etbin],"L0RMEQ", "", fitResoEbinsEta[i][ebin][etbin]->GetParameter(1)-1*fitResoEbinsEta[i][ebin][etbin]->GetParameter(2), fitResoEbinsEta[i][ebin][etbin]->GetParameter(1)+1.5*fitResoEbinsEta[i][ebin][etbin]->GetParameter(2));
//                 cout << "particle: " << i <<  " eta bin: " << etbin << " \t ebin: "<< ebin<< "\t" << fitResoEbinsEta[i][ebin][etbin]->GetChisquare()/fitResoEbinsEta[i][ebin][etbin]->GetNDF() << endl;
                meanEta     = fitResoEbinsEta[i][ebin][etbin]->GetParameter(1);
                meanEtaErr  = fitResoEbinsEta[i][ebin][etbin]->GetParError(1);
                sigmaEta    = fitResoEbinsEta[i][ebin][etbin]->GetParameter(2);
                sigmaEtaErr = fitResoEbinsEta[i][ebin][etbin]->GetParError(2);
                if (i > 1){
                  fwhmEta     = sigmaEta;
                  fwhmEtaErr  = sigmaEtaErr;
                } else {
                  fwhmEta  = (fitResoEbinsEta[i][ebin][etbin]->GetParameter(1) - fitResoEbinsEta[i][ebin][etbin]->GetX(0.5*fitResoEbinsEta[i][ebin][etbin]->GetMaximum(), 0.2 , fitResoEbinsEta[i][ebin][etbin]->GetParameter(1)));
                  Double_t relEM = fitResoEbinsEta[i][ebin][etbin]->GetParError(1)/ fitResoEbinsEta[i][ebin][etbin]->GetParameter(1);
                  Double_t relES = fitResoEbinsEta[i][ebin][etbin]->GetParError(2)/ fitResoEbinsEta[i][ebin][etbin]->GetParameter(2);
                  Double_t relEL = fitResoEbinsEta[i][ebin][etbin]->GetParError(3)/ fitResoEbinsEta[i][ebin][etbin]->GetParameter(3);
                  Double_t relEN = fitResoEbinsEta[i][ebin][etbin]->GetParError(4)/ fitResoEbinsEta[i][ebin][etbin]->GetParameter(4);
                  fwhmEtaErr = fwhmEta*TMath::Sqrt(relEM*relEM+relES*relES+relEL*relEL+relEN*relEN);
                }
              }
              Double_t resoEta     = sigmaEta/meanEta;
              Double_t resoEtaErr  = TMath::Sqrt(TMath::Power(sigmaEtaErr/sigmaEta,2)+TMath::Power(meanEtaErr/meanEta,2));
              Double_t resoFEta     = fwhmEta/meanEta;
              Double_t resoFEtaErr  = TMath::Sqrt(TMath::Power(fwhmEtaErr/fwhmEta,2)+TMath::Power(meanEtaErr/meanEta,2));

              histResoEVsEEta[i][etbin]->SetBinContent(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), 100*resoEta);
              histResoEVsEEta[i][etbin]->SetBinError(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), resoEtaErr);
              histResoEVs1oEEta[i][etbin]->SetBinContent(histResoEVs1oEEta[i][etbin]->FindBin(1./TMath::Sqrt(partE[ebin])), 100*resoEta);
              histResoEVs1oEEta[i][etbin]->SetBinError(histResoEVs1oEEta[i][etbin]->FindBin(1./TMath::Sqrt(partE[ebin])), resoEtaErr);
              histResoFEVsEEta[i][etbin]->SetBinContent(histResoFEVsEEta[i][etbin]->FindBin(partE[ebin]), 100*resoFEta);
              histResoFEVsEEta[i][etbin]->SetBinError(histResoFEVsEEta[i][etbin]->FindBin(partE[ebin]), resoFEtaErr);
              histResoFEVs1oEEta[i][etbin]->SetBinContent(histResoFEVs1oEEta[i][etbin]->FindBin(1./TMath::Sqrt(partE[ebin])), 100*resoFEta);
              histResoFEVs1oEEta[i][etbin]->SetBinError(histResoFEVs1oEEta[i][etbin]->FindBin(1./TMath::Sqrt(partE[ebin])), resoFEtaErr);

              histMeanEVsEEta[i][etbin]->SetBinContent(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), meanEta);
              histMeanEVsEEta[i][etbin]->SetBinError(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), meanEtaErr);
              histSigmaEVsEEta[i][etbin]->SetBinContent(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), sigmaEta);
              histSigmaEVsEEta[i][etbin]->SetBinError(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), sigmaEtaErr);
              histFWHMEVsEEta[i][etbin]->SetBinContent(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), fwhmEta);
              histFWHMEVsEEta[i][etbin]->SetBinError(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), fwhmEtaErr);
            } else {
              useEEta[i][ebin][etbin] = 0;
            }
          }
        } else {
          for (Int_t ebin = 0; ebin < nEne; ebin++){
            useEEta[i][ebin][etbin] = 0;
          } 
        }
      }
    } else {
      for (Int_t ebin = 0; ebin < nEne; ebin++){
        useE[i][ebin] = 0;
      }
    } 
  }
  
  cout << "available energies & particles" << endl;
  cout << "E\t" ;
  for (Int_t i = 0; i < nParticles; i++){
     cout <<  labelPart[i] << "\t";
     minEPart[i] = 1000;
     maxEPart[i] = 0;
  }
  cout << endl;
  for (Int_t ebin = 0; ebin < nEne; ebin++){
    cout << partE[ebin] << "\t" ;
    for (Int_t i = 0; i < nParticles; i++){
      cout <<  useE[i][ebin] << "\t";
      if (useE[i][ebin]){
        if (minEPart[i] > partE[ebin]) minEPart[i] = partE[ebin];
        if (maxEPart[i] < partE[ebin]) maxEPart[i] = partE[ebin];
      }
    }
    cout << endl;
  }

  for (Int_t i = 0; i < nParticles; i++){
    if (usePart[i]){
      nActivePart++;
      cout <<  labelPart[i] << "\t min E \t" << minEPart[i] << "\t max E\t"<< maxEPart[i] << endl ;
    }
  }
  TLegend* legendResYRECAL    = GetAndSetLegend2(0.14, 0.14, 0.6, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYRECAL->AddEntry(graphReqEM, "YR requirement EM","f");
  legendResYRECAL->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermECal[dirCal][0], sqrtETermECal[dirCal][1], linTermECal[dirCal][0], linTermECal[dirCal][1]),"");
  TLegend* legendResYR1oEECal = GetAndSetLegend2(0.58, 0.14, 0.97, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYR1oEECal->AddEntry(graphReq1oEEM, "YR requirement EM","f");
  legendResYR1oEECal->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermECal[dirCal][0], sqrtETermECal[dirCal][1], linTermECal[dirCal][0], linTermECal[dirCal][1]),"");
  TLegend* legendResYRHCAL    = GetAndSetLegend2(0.14, 0.14, 0.6, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYRHCAL->AddEntry(graphReqHad, "YR requirement Had","f");
  legendResYRHCAL->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermHCal[dirCal][0], sqrtETermHCal[dirCal][1], linTermHCal[dirCal][0], linTermHCal[dirCal][1]),"");
  TLegend* legendResYR1oEHCal = GetAndSetLegend2(0.58, 0.14, 0.97, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYR1oEHCal->AddEntry(graphReq1oEHad, "YR requirement Had","f");
  legendResYR1oEHCal->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermHCal[dirCal][0], sqrtETermHCal[dirCal][1], linTermHCal[dirCal][0], linTermHCal[dirCal][1]),"");

  // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  for (Int_t i = 0; i < nParticles-1; i++){
    if (!usePart[i]) continue;
    TCanvas* cPNG;
    split_canvas(cPNG, Form("cPNG%d",i), nActiveEpart[i]);
    Int_t padnum =0;
    std::cout << labelPart[i] << std::endl;
    Float_t maxReso = maxResHad[dirCal];
    if (i < 2) 
      maxReso = maxResEM[dirCal];
    for(Int_t ebin=0; ebin< nEne; ebin++){
      if (!useE[i][ebin]) continue;
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
      cPNG->cd(padnum+1)->SetLogy(1);
      SetStyleHistoTH1ForGraphs(hist1DResoEbins[i][ebin], "#it{E}^{rec}/#it{E}^{MC}", "counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
      hist1DResoEbins[i][ebin]->GetXaxis()->SetRangeUser(0.,maxReso);
      hist1DResoEbins[i][ebin]->GetYaxis()->SetRangeUser(0.8,hist1DResoEbins[i][ebin]->GetMaximum()*2);
      
      DrawGammaSetMarker(hist1DResoEbins[i][ebin], 20, 1.5, kBlack, kBlack);      
      hist1DResoEbins[i][ebin]->Draw("p,e");
      DrawGammaSetMarker(hist1DResoEbinshighests[i][ebin], 24, 1.5, kGreen+2, kGreen+2);      
      hist1DResoEbinshighests[i][ebin]->Draw("p,e,same");
      
      if (fitResoEbins[i][ebin] && doFitting){
        DrawGammaSetMarkerTF1( fitResoEbins[i][ebin], 7, 2, kBlue+1);
        fitResoEbins[i][ebin]->Draw("same");
      }
      if (fitResoEbinshighest[i][ebin] && doFitting){
        DrawGammaSetMarkerTF1( fitResoEbinshighest[i][ebin], 9, 2, kGreen+3);
        fitResoEbinshighest[i][ebin]->Draw("same");
      }
      
      if (padnum == 0){
        drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.1f< #eta< %1.1f",minEtaMax,maxEtaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
       
      drawLatexAdd(Form("#it{E} (GeV) = %1.1f ",partE[ebin]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      Int_t bin = histMeanEVsE[i]->FindBin(partE[ebin]);
      DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin), 
                     histMeanEVsE[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoEbins[i][ebin]->GetMaximum(), 1, kGray+2, 2);
      if (i < 2){
        DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin)-histFWHMEVsE[i]->GetBinContent(bin), 
                      histMeanEVsE[i]->GetBinContent(bin)-histFWHMEVsE[i]->GetBinContent(bin), 
                      0., 0.05*hist1DResoEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
      } else {
         DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin)-histSigmaEVsE[i]->GetBinContent(bin), 
                      histMeanEVsE[i]->GetBinContent(bin)-histSigmaEVsE[i]->GetBinContent(bin), 
                      0., 0.05*hist1DResoEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
      }
      DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin)+histSigmaEVsE[i]->GetBinContent(bin), 
                    histMeanEVsE[i]->GetBinContent(bin)+histSigmaEVsE[i]->GetBinContent(bin), 
                    0., 0.05*hist1DResoEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
      padnum++;
    }
    cPNG->Print(Form("%s/Energy_%sResolution.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));

    for (Int_t etbin = 0; etbin < nEta; etbin++){
      if (!useEta[i][etbin]) continue;
      Double_t etaMinCurr = partEtaCalo[etbin];
      Double_t etaMaxCurr = partEtaCalo[etbin+1];
      split_canvas(cPNG, Form("cPNG%d_%d",i,etbin), nActiveEpart[i]);
      Int_t padnum =0;
      std::cout << etaMinCurr << "\t" << etaMaxCurr << std::endl;
        
      for(Int_t ebin=0; ebin< nEne; ebin++){
        if (!useEEta[i][ebin][etbin]) continue;
        cPNG->cd(padnum+1);
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1)->SetLogy(1);
        SetStyleHistoTH1ForGraphs(hist1DResoEbinsEta[i][ebin][etbin], "#it{E}^{rec}/#it{E}^{MC}", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        hist1DResoEbinsEta[i][ebin][etbin]->GetXaxis()->SetRangeUser(0.,maxReso);
        hist1DResoEbinsEta[i][ebin][etbin]->GetYaxis()->SetRangeUser(0.8,hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum()*2);
        
        DrawGammaSetMarker(hist1DResoEbinsEta[i][ebin][etbin], 20, 1.5, kBlack, kBlack);      
        hist1DResoEbinsEta[i][ebin][etbin]->Draw("p,e");
        
        if (fitResoEbinsEta[i][ebin][etbin] && doFitting ){
          DrawGammaSetMarkerTF1( fitResoEbinsEta[i][ebin][etbin], 7, 2, kBlue+1);
          fitResoEbinsEta[i][ebin][etbin]->Draw("same");
        }
        
        if (padnum == 0){
          drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%1.1f< #eta< %1.1f",etaMinCurr,etaMaxCurr),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        
        drawLatexAdd(Form("#it{E} (GeV) = %1.1f ",partE[ebin]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        Int_t bin = histMeanEVsEEta[i][etbin]->FindBin(partE[ebin]);
        DrawGammaLines(histMeanEVsEEta[i][etbin]->GetBinContent(bin), 
                      histMeanEVsEEta[i][etbin]->GetBinContent(bin), 
                      0., 0.05*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 1, kGray+2, 2);
        if (i < 2){
          DrawGammaLines(histMeanEVsEEta[i][etbin]->GetBinContent(bin)-histSigmaEVsEEta[i][etbin]->GetBinContent(bin), 
                        histMeanEVsEEta[i][etbin]->GetBinContent(bin)-histSigmaEVsEEta[i][etbin]->GetBinContent(bin), 
                        0., 0.05*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 1, kRed+2, 2);
        }else {
          DrawGammaLines(histMeanEVsEEta[i][etbin]->GetBinContent(bin)-histFWHMEVsEEta[i][etbin]->GetBinContent(bin), 
                        histMeanEVsEEta[i][etbin]->GetBinContent(bin)-histFWHMEVsEEta[i][etbin]->GetBinContent(bin), 
                        0., 0.05*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 1, kRed+2, 2);
          
        }
        DrawGammaLines(histMeanEVsEEta[i][etbin]->GetBinContent(bin)+histSigmaEVsEEta[i][etbin]->GetBinContent(bin), 
                      histMeanEVsEEta[i][etbin]->GetBinContent(bin)+histSigmaEVsEEta[i][etbin]->GetBinContent(bin), 
                      0., 0.05*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNG->Print(Form("%s/Energy_%sResolution_Etabin%d.%s", outputDir.Data(), readNames[i].Data(), etbin, suffix.Data()));

    }
    
    if (resFit1oEEM[dirCal][1] ==  -1) resFit1oEEM[dirCal][1]  = 1/TMath::Sqrt(minEPart[i]);
    if (resFit1oEEM[dirCal][0] ==  -1) resFit1oEEM[dirCal][0]  = 1/TMath::Sqrt(maxEPart[i]);
    if (resFitEM[dirCal][1] ==  -1) resFitEM[dirCal][1]    = maxEPart[i];
    if (resFitEM[dirCal][0] ==  -1) resFitEM[dirCal][0]    = minEPart[i];
    if (meanFitEM[dirCal][1] ==  -1) meanFitEM[dirCal][1]  = maxEPart[i];
    if (meanFitEM[dirCal][0] ==  -1) meanFitEM[dirCal][0]  = minEPart[i];

    //*************************************************************************************
    // WIDTH PLOT
    //*************************************************************************************
    TCanvas* cWidth = new TCanvas("cWidth","",0,0,1100,1000);
    cWidth->SetLogx(1);
    DrawGammaCanvasSettings( cWidth, 0.095, 0.01, 0.01, 0.105);
    TH2F * histWidthDummy                           = new TH2F("histWidthDummy","histWidthDummy",1000,0.3, 170,1000,0, maxResHad[dirCal]-1);
    SetStyleHistoTH2ForGraphs(histWidthDummy, "#it{E} (GeV)","#sigma (GeV)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    histWidthDummy->GetXaxis()->SetNoExponent();
    histWidthDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histWidthDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    if (i < 2) 
      histWidthDummy->GetYaxis()->SetRangeUser(0, maxResEM[dirCal]-1);
    else 
      histWidthDummy->GetYaxis()->SetRangeUser(0, maxResHad[dirCal]-1);
    histWidthDummy->Draw();

    DrawGammaSetMarker(histSigmaEVsE[i], markerStylePart[i], markerSizePart[i], colorPart[i], colorPart[i]);
    histSigmaEVsE[i]->Draw("same,p");
    if (i < 2){
      DrawGammaSetMarker(histFWHMEVsE[i], markerStylePart[i], markerSizePart[i], colorPart[i]-6, colorPart[i]-6);
      histFWHMEVsE[i]->Draw("same,p");
    }
    
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    cWidth->Print(Form("%s/Width_%sFitted.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));

    DrawGammaSetMarker(histSigmaEVsEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]+2, colorPart[i]+2);
    histSigmaEVsEhighest[i]->Draw("same,p");
        if (i < 2){
      DrawGammaSetMarker(histFWHMEVsEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]-4, colorPart[i]-4);
      histFWHMEVsEhighest[i]->Draw("same,p");
    }

    cWidth->Print(Form("%s/Width_%sFittedWithHighest.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));

    //*************************************************************************************
    // RESOLUTION PLOT
    //*************************************************************************************
    TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
    cReso->SetLogx(1);
    DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
    TH2F * histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0.3, 170,1000,0, maxResPerHad[dirCal]);
    SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    histResoDummy->GetXaxis()->SetNoExponent();
    histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    if (i < 2) 
      histResoDummy->GetYaxis()->SetRangeUser(0, maxResPerEM[dirCal]);
    else 
      histResoDummy->GetYaxis()->SetRangeUser(0, maxResPerHad[dirCal]);
    histResoDummy->Draw();

    TLegend* legendResLabE  = GetAndSetLegend2(0.14, 0.95-0.85*textSizeLabelsRel, 0.35, 0.95,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);

    DrawGammaSetMarker(histResoEVsE[i], markerStylePart[i], markerSizePart[i], colorPart[i], colorPart[i]);
    histResoEVsE[i]->Draw("same,p");
//     fitResoE[i] = new TF1(Form("fit_reso_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x)))", minEPart[i],maxEPart[i]);
    fitResoE[i] = new TF1(Form("fit_reso_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x))+([2]/x)", resFitEM[dirCal][0], resFitEM[dirCal][1]);
    fitResoE[i]->SetParLimits(0,0.,20);
    fitResoE[i]->SetParLimits(1,0.,90);
    fitResoE[i]->SetParameter(2,0.01);
    fitResoE[i]->SetParLimits(2,0.,0.3);
    histResoEVsE[i]->Fit(fitResoE[i],"RMNE");
    fitResoE[i]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fitResoE[i], 3, 3, colorPart[i]+1);
    fitResoE[i]->SetRange(minEPart[i], maxEPart[i]);
    fitResoE[i]->Draw("same");
    if (i < 2){
      DrawGammaSetMarker(histResoFEVsE[i], markerStylePart[i], markerSizePart[i], colorPart[i]-6, colorPart[i]-6);
      histResoFEVsE[i]->Draw("same,p");
      fitResoFE[i] = new TF1(Form("fit_resof_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x))+([2]/x)", resFitEM[dirCal][0], resFitEM[dirCal][1]);
      fitResoFE[i]->SetParLimits(0,0.,20);
      fitResoFE[i]->SetParLimits(1,0.,90);
      fitResoFE[i]->SetParameter(2,0.01);
      fitResoFE[i]->SetParLimits(2,0.,0.3);
      histResoFEVsE[i]->Fit(fitResoFE[i],"RMNE");
      fitResoFE[i]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fitResoFE[i], 3, 3, colorPart[i]-6);
      fitResoFE[i]->SetRange(minEPart[i], maxEPart[i]);
      fitResoFE[i]->Draw("same");
      
    }
    
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendResLabE->AddEntry(fitResoE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
    legendResLabE->Draw();

    cReso->Print(Form("%s/Resolution_%sFitted.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));

    DrawGammaSetMarker(histResoEVsEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]+2, colorPart[i]+2);
    histResoEVsEhighest[i]->Draw("same,p");
    fitResoEhighest[i] = new TF1(Form("fit_resohighest_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x))+([2]/x)", resFitEM[dirCal][0],resFitEM[dirCal][1]);
    fitResoEhighest[i]->SetParLimits(0,0.,20);
    fitResoEhighest[i]->SetParLimits(1,0.,90);
    fitResoEhighest[i]->SetParameter(2,0.01);
    fitResoEhighest[i]->SetParLimits(2,0.,0.3);
    histResoEVsEhighest[i]->Fit(fitResoEhighest[i],"RMNE");
    DrawGammaSetMarkerTF1(fitResoEhighest[i], 7, 3, colorPart[i]+2);
    fitResoEhighest[i]->SetRange(minEPart[i], maxEPart[i]);
    fitResoEhighest[i]->Draw("same");
    if (i < 2){
      DrawGammaSetMarker(histResoFEVsEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]-4, colorPart[i]-4);
      histResoFEVsEhighest[i]->Draw("same,p");
      fitResoFEhighest[i] = new TF1(Form("fit_resohighestf_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x))+([2]/x)", resFitEM[dirCal][0], resFitEM[dirCal][1]);
      fitResoFEhighest[i]->SetParLimits(0,0.,20);
      fitResoFEhighest[i]->SetParLimits(1,0.,90);
      fitResoFEhighest[i]->SetParameter(2,0.01);
      fitResoFEhighest[i]->SetParLimits(2,0.,0.3);
      histResoFEVsEhighest[i]->Fit(fitResoFEhighest[i],"RMNE");
      fitResoFEhighest[i]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fitResoFEhighest[i], 3, 3, colorPart[i]-4);
      fitResoFEhighest[i]->SetRange(minEPart[i], maxEPart[i]);
      fitResoFEhighest[i]->Draw("same");
    }

    legendResLabE  = GetAndSetLegend2(0.14, 0.95-2*0.85*textSizeLabelsRel, 0.35, 0.95,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);
    legendResLabE->AddEntry(fitResoE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
    legendResLabE->AddEntry(fitResoEhighest[i], Form("%s highest: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitResoEhighest[i]->GetParameter(1),fitResoEhighest[i]->GetParameter(0)),"l");
    legendResLabE->Draw();

    
    cReso->Print(Form("%s/Resolution_%sFittedWithHighest.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
    histResoDummy->Draw();
    TLegend* legendResEEta  = GetAndSetLegend2(0.14, 0.95-((nActiveEtapart[i]+1)*0.85*textSizeLabelsRel), 0.35, 0.95,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);

    if (i < 2)
      graphReqEM->Draw("same,e3");
    else 
      graphReqHad->Draw("same,e3");
    
    histResoDummy->Draw("axis,same");
    for (Int_t etbin = minEtaBinCaloDis[dirCal]; etbin < maxEtaBinCaloDis[dirCal]; etbin++){
      if (!useEta[i][etbin]) continue;
      DrawGammaSetMarker(histResoEVsEEta[i][etbin], markerStyleEta[etbin], markerSizeEta[etbin]*1.5, colorEta[etbin], colorEta[etbin]);
      histResoEVsEEta[i][etbin]->Draw("same,p");
      legendResEEta->AddEntry(histResoEVsEEta[i][etbin],Form("%1.1f< #eta< %1.1f", partEtaCalo[etbin], partEtaCalo[etbin+1]),"p");
    }
    DrawGammaSetMarkerTF1(fitResoE[i], 3, 3, kBlack);
    fitResoE[i]->Draw("same");
    legendResEEta->AddEntry(fitResoE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
    if (i < 2)
      legendResYRECAL->Draw();
    else 
      legendResYRHCAL->Draw();
    
    legendResEEta->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/Resolution_%sFitted_etaDep.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
      
  //*************************************************************************************
  // RESOLUTION PLOT 2
  //*************************************************************************************

    TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
    DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
    histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 1/TMath::Sqrt(minEPart[i])+0.1,1000,0, maxResPerHad[dirCal]);
    SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    histResoDummy->GetXaxis()->SetNoExponent();
    histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    histResoDummy->Draw();
    
    if (i < 2) 
      histResoDummy->GetYaxis()->SetRangeUser(0, maxResPerEM[dirCal]);
    else 
      histResoDummy->GetYaxis()->SetRangeUser(0, maxResPerHad[dirCal]);

    
    TLegend* legendResLab1oE  = GetAndSetLegend2(0.14, 0.91-0.85*textSizeLabelsRel, 0.35, 0.91,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);
    DrawGammaSetMarker(histResoEVs1oE[i], markerStylePart[i], markerSizePart[i], colorPart[i], colorPart[i]);
    histResoEVs1oE[i]->Draw("same,p");
    fitReso1oE[i] = new TF1(Form("fit_reso_%s1oE", readNames[i].Data()), "[0]+[1]*x+[2]*TMath::Power(x,2)", resFit1oEEM[dirCal][0],resFit1oEEM[dirCal][1]); //+[2]*TMath::Power(x,2)
    fitReso1oE[i]->SetParameter(2,0.01);
    fitReso1oE[i]->SetParLimits(2,0.001, 0.3);

    histResoEVs1oE[i]->Fit(fitReso1oE[i],"RMNE");
    fitReso1oE[i]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fitReso1oE[i], 3, 3, colorPart[i]+1);
    fitReso1oE[i]->SetRange(1/TMath::Sqrt(maxEPart[i]), 1/TMath::Sqrt(minEPart[i]));
    fitReso1oE[i]->Draw("same");
    
    if (i < 2){
      DrawGammaSetMarker(histResoFEVs1oE[i], markerStylePart[i], markerSizePart[i], colorPart[i]-6, colorPart[i]-6);
      histResoFEVs1oE[i]->Draw("same,p");
      fitResoF1oE[i] = new TF1(Form("fit_resof_%s1oE", readNames[i].Data()), "[0]+[1]*x+[2]*TMath::Power(x,2)", resFit1oEEM[dirCal][0],resFit1oEEM[dirCal][1]); //+[2]*TMath::Power(x,2)
      fitResoF1oE[i]->SetParameter(2,0.01);
      fitResoF1oE[i]->SetParLimits(2,0.001, 0.3);

      histResoFEVs1oE[i]->Fit(fitResoF1oE[i],"RMNE");
      fitResoF1oE[i]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fitResoF1oE[i], 3, 3, colorPart[i]-6);
      fitResoF1oE[i]->SetRange(1/TMath::Sqrt(maxEPart[i]), 1/TMath::Sqrt(minEPart[i]));
      fitResoF1oE[i]->Draw("same"); 
    }

    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.14,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    legendResLab1oE->AddEntry(fitReso1oE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitReso1oE[i]->GetParameter(1),fitReso1oE[i]->GetParameter(0)),"l");
    legendResLab1oE->Draw();

    cReso2->Print(Form("%s/Resolution_%sFitted_1oE.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));


    DrawGammaSetMarker(histResoEVs1oEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]+2, colorPart[i]+2);
    histResoEVs1oEhighest[i]->Draw("same,p");
    fitReso1oEhighest[i] = new TF1(Form("fit_resohighest_%s1oE", readNames[i].Data()), "[0]+[1]*x+[2]*TMath::Power(x,2)", resFit1oEEM[dirCal][0],resFit1oEEM[dirCal][1]); //+[2]*TMath::Power(x,2)
    fitReso1oEhighest[i]->SetParameter(2,0.01);
    fitReso1oEhighest[i]->SetParLimits(2,0.0, 0.3);
    histResoEVs1oEhighest[i]->Fit(fitReso1oEhighest[i],"RMNE");
    DrawGammaSetMarkerTF1(fitReso1oEhighest[i], 7, 3, colorPart[i]+2);
    fitReso1oEhighest[i]->SetRange(1/TMath::Sqrt(maxEPart[i]), 1/TMath::Sqrt(minEPart[i]));
    fitReso1oEhighest[i]->Draw("same");
    
    if (i < 2){
      DrawGammaSetMarker(histResoFEVs1oEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]-4, colorPart[i]-4);
      histResoFEVs1oEhighest[i]->Draw("same,p");
      fitResoF1oEhighest[i] = new TF1(Form("fit_resofhighest_%s1oE", readNames[i].Data()), "[0]+[1]*x+[2]*TMath::Power(x,2)", resFit1oEEM[dirCal][0],resFit1oEEM[dirCal][1]); //+[2]*TMath::Power(x,2)
      fitResoF1oEhighest[i]->SetParameter(2,0.01);
      fitResoF1oEhighest[i]->SetParLimits(2,0.0, 0.3);
      histResoFEVs1oEhighest[i]->Fit(fitResoF1oEhighest[i],"RMNE");
      DrawGammaSetMarkerTF1(fitResoF1oEhighest[i], 7, 3, colorPart[i]-4);
      fitResoF1oEhighest[i]->SetRange(1/TMath::Sqrt(maxEPart[i]), 1/TMath::Sqrt(minEPart[i]));
      fitResoF1oEhighest[i]->Draw("same");
    
    }  
      
    legendResLab1oE  = GetAndSetLegend2(0.14, 0.91-2*0.85*textSizeLabelsRel, 0.35, 0.91,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);
    legendResLab1oE->AddEntry(fitReso1oE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitReso1oE[i]->GetParameter(1),fitReso1oE[i]->GetParameter(0)),"l");
    legendResLab1oE->AddEntry(fitReso1oEhighest[i], Form("%s highest: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitReso1oEhighest[i]->GetParameter(1),fitReso1oEhighest[i]->GetParameter(0)),"l");
    legendResLab1oE->Draw();

    
    cReso2->Print(Form("%s/Resolution_%sFittedWithHighest_1oE.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
    
    histResoDummy->Draw();
    if (i < 2)
      graphReq1oEEM->Draw("same,e3");
    else 
      graphReq1oEHad->Draw("same,e3");
    
    histResoDummy->Draw("axis,same");
    for (Int_t etbin = minEtaBinCaloDis[dirCal]; etbin < maxEtaBinCaloDis[dirCal]; etbin++){
      if (!useEta[i][etbin]) continue;
      DrawGammaSetMarker(histResoEVs1oEEta[i][etbin], markerStyleEta[etbin], markerSizeEta[etbin]*1.5, colorEta[etbin], colorEta[etbin]);
      histResoEVs1oEEta[i][etbin]->Draw("same,p");
    }
    DrawGammaSetMarkerTF1(fitReso1oE[i], 3, 3, kBlack);
    fitReso1oE[i]->Draw("same");
    if (i < 2)
      legendResYR1oEECal->Draw();
    else 
      legendResYR1oEHCal->Draw();
    
    legendResEEta->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    cReso2->Print(Form("%s/Resolution_%sFitted_1oE_etaDep.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
  }

  if (isEMCal && nActivePart >2 ) nActivePart = 2;
    else nActivePart = nActivePart-1;

  cout << nActivePart << endl;
  TLegend* legendResYRCAL    = GetAndSetLegend2(0.56, 0.86-(nActivePart*0.85*textSizeLabelsRel), 0.97, 0.86-((4+nActivePart)*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYRCAL->AddEntry(graphReqEM, "YR requirement EM","f");
  legendResYRCAL->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermECal[dirCal][0], sqrtETermECal[dirCal][1], linTermECal[dirCal]    
  [0], linTermECal[dirCal][1]),"");
  legendResYRCAL->AddEntry(graphReqHad, "YR requirement Had","f");
  legendResYRCAL->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermHCal[dirCal][0], sqrtETermHCal[dirCal][1], linTermHCal[dirCal][0], linTermECal[dirCal][1]),"");
  TLegend* legendResYR1oECAL = GetAndSetLegend2(0.14,0.95-((nActivePart+1)*0.85*textSizeLabelsRel), 0.6, 0.95-((4+nActivePart+1)*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYR1oECAL->AddEntry(graphReq1oEEM, "YR requirement EM","f");
  legendResYR1oECAL->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermECal[dirCal][0], sqrtETermECal[dirCal][1], linTermECal[dirCal][0], linTermECal[dirCal][1]),"");
  legendResYR1oECAL->AddEntry(graphReq1oEHad, "YR requirement Had","f");
  legendResYR1oECAL->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETermHCal[dirCal][0], sqrtETermHCal[dirCal][1], linTermHCal[dirCal][0], linTermHCal[dirCal][1]),"");

  
  TCanvas* cExampleBin = new TCanvas("cExampleBin","",0,0,1100,1000);
  DrawGammaCanvasSettings( cExampleBin, 0.095, 0.01, 0.045, 0.105);
  
    Float_t maximYRange = 0;
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!useE[i][exEbin]) continue;
      if (maximYRange < hist1DResoEbinshighests[i][exEbin]->GetMaximum())
        maximYRange = hist1DResoEbinshighests[i][exEbin]->GetMaximum();
    }
  
  
    TLegend* legendExamBin  = GetAndSetLegend2(0.14, 0.77-((nActivePart+2)/2*0.85*textSizeLabelsRel), 0.35, 0.77,0.85*textSizeLabelsPixel, 2, "", 43, 0.35);
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!useE[i][exEbin]) continue;
      SetStyleHistoTH1ForGraphs(hist1DResoEbinshighests[i][exEbin], "#it{E}^{rec}/#it{E}^{MC}", "counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
      hist1DResoEbinshighests[i][exEbin]->GetXaxis()->SetRangeUser(0.,maxReso);
      hist1DResoEbinshighests[i][exEbin]->GetYaxis()->SetRangeUser(0.,1.2*maximYRange);
      
      DrawGammaSetMarker(hist1DResoEbinshighests[i][exEbin], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      if (i == 0)
        hist1DResoEbinshighests[i][exEbin]->Draw("p,e");
      else 
        hist1DResoEbinshighests[i][exEbin]->Draw("p,e,same");
      legendExamBin->AddEntry(hist1DResoEbinshighests[i][exEbin], labelPart[i], "p");
      
      if (fitResoEbinshighest[i][exEbin] && doFitting ){
        DrawGammaSetMarkerTF1( fitResoEbinshighest[i][exEbin], lineStylePart[i], 2, colorPart[i]);
        fitResoEbinshighest[i][exEbin]->Draw("same");
      }
    }
    legendExamBin->Draw();
    drawLatexAdd(perfLabel,0.14,0.88,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.14,0.83,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("%1.1f< #eta< %1.1f, #it{E} = %1.1f GeV", minEtaMax, maxEtaMax, partE[exEbin]),0.14,0.78,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    
  cExampleBin->Print(Form("%s/ExampleBin.%s", outputDir.Data(), suffix.Data()));
  
  
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
  cReso->SetLogx(1);
  DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
  TH2F * histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0.3, 170,1000,0, maxResPerHad[dirCal]);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);  
  if (isEMCal) 
    histResoDummy->GetYaxis()->SetRangeUser(0,maxResPerEM[dirCal]);
  else 
    histResoDummy->GetYaxis()->SetRangeUser(0,maxResPerHad[dirCal]);
  histResoDummy->Draw();
  
    if (doFitting && usePart[0])  nActivePart++;
    if (doFitting && usePart[1])  nActivePart++;
    TLegend* legendResPart  = GetAndSetLegend2(0.58, 0.86-(nActivePart*0.85*textSizeLabelsRel), 0.97, 0.86,0.85*textSizeLabelsPixel, 2, "", 43, 0.22);
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      DrawGammaSetMarker(histResoEVsE[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histResoEVsE[i]->Draw("same,p");
      DrawGammaSetMarkerTF1(fitResoE[i], lineStylePart[i], 3, colorPart[i]+1);
      fitResoE[i]->Draw("same");
      if (i < 2 && doFitting){
        DrawGammaSetMarker(histResoFEVsE[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i]-6, colorPart[i]-6);
        histResoFEVsE[i]->Draw("same,p");
        DrawGammaSetMarkerTF1(fitResoFE[i], lineStylePart[i], 3, colorPart[i]-6);
        fitResoFE[i]->Draw("same");
        legendResPart->AddEntry(histResoEVsE[i],Form("%s, r", labelPart[i].Data()),"p");
        legendResPart->AddEntry(fitResoE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
        legendResPart->AddEntry(histResoFEVsE[i],Form("%s, l", labelPart[i].Data()),"p");
        legendResPart->AddEntry(fitResoFE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoFE[i]->GetParameter(1),fitResoFE[i]->GetParameter(0)),"l");
      } else {      
        legendResPart->AddEntry(histResoEVsE[i],Form("%s", labelPart[i].Data()),"p");
        legendResPart->AddEntry(fitResoE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
      }

    }
    legendResPart->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/Resolution_Fitted_allPart.%s", outputDir.Data(), suffix.Data()));

  histResoDummy->Draw();
    TLegend* legendResPartHighest  = GetAndSetLegend2(0.58, 0.86-(nActivePart*0.85*textSizeLabelsRel), 0.97, 0.86,0.85*textSizeLabelsPixel, 2, "", 43, 0.22); 
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      DrawGammaSetMarker(histResoEVsEhighest[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histResoEVsEhighest[i]->Draw("same,p");
      DrawGammaSetMarkerTF1(fitResoEhighest[i], lineStylePart[i], 3, colorPart[i]+1);
      fitResoEhighest[i]->Draw("same");
      
      if (i < 2 && doFitting){
        DrawGammaSetMarker(histResoFEVsEhighest[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i]-6, colorPart[i]-6);
        histResoFEVsEhighest[i]->Draw("same,p");
        DrawGammaSetMarkerTF1(fitResoFEhighest[i], lineStylePart[i], 3, colorPart[i]-6);
        fitResoFEhighest[i]->Draw("same");
        legendResPartHighest->AddEntry(histResoEVsEhighest[i],Form("%s, r", labelPart[i].Data()),"p");
        legendResPartHighest->AddEntry(fitResoEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoEhighest[i]->GetParameter(1),fitResoEhighest[i]->GetParameter(0)),"l");
        legendResPartHighest->AddEntry(histResoFEVsEhighest[i],Form("%s, l", labelPart[i].Data()),"p");
        legendResPartHighest->AddEntry(fitResoFEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoFEhighest[i]->GetParameter(1),fitResoFEhighest[i]->GetParameter(0)),"l");
      } else {      
        legendResPartHighest->AddEntry(histResoEVsEhighest[i],Form("%s", labelPart[i].Data()),"p");
        legendResPartHighest->AddEntry(fitResoEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoEhighest[i]->GetParameter(1),fitResoEhighest[i]->GetParameter(0)),"l");
      }
    }
    legendResPartHighest->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/Resolution_FittedHighest_allPart.%s", outputDir.Data(), suffix.Data()));

  histResoDummy->Draw();
    if (isEMCal){
      graphReqEM->Draw("same,e3");
    } else { 
      graphReqEM->Draw("same,e3");
      graphReqHad->Draw("same,e3");
    }
      
    histResoDummy->Draw("axis,same");
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      histResoEVsEhighest[i]->Draw("same,p");
      fitResoEhighest[i]->Draw("same");
      if (i < 2 && doFitting){
        histResoFEVsEhighest[i]->Draw("same,p");
        fitResoFEhighest[i]->Draw("same");
      }
      
    }
    legendResPartHighest->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (isEMCal)
      legendResYRECAL->Draw();
    else 
      legendResYRCAL->Draw();
    
  cReso->Print(Form("%s/Resolution_FittedHighest_allPart_wReq.%s", outputDir.Data(), suffix.Data()));

  
// RESOLUTION PLOT 2

  TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 1.5,1000,0, maxResPerHad[dirCal]);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    if (isEMCal) 
    histResoDummy->GetYaxis()->SetRangeUser(0,maxResPerEM[dirCal]);
  else 
    histResoDummy->GetYaxis()->SetRangeUser(0,maxResPerHad[dirCal]);
  histResoDummy->Draw();
  
  
  
    TLegend* legendResPart1oE  = GetAndSetLegend2(0.14, 0.95-(nActivePart*0.85*textSizeLabelsRel), 0.6, 0.95,0.85*textSizeLabelsPixel, 2, "", 43, 0.2);
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      DrawGammaSetMarker(histResoEVs1oE[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histResoEVs1oE[i]->Draw("same,p");
      DrawGammaSetMarkerTF1(fitReso1oE[i], lineStylePart[i], 3, colorPart[i]+1);
      fitReso1oE[i]->Draw("same");
      
      if (i < 2 && doFitting){
        DrawGammaSetMarker(histResoFEVs1oE[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i]-6, colorPart[i]-6);
        histResoFEVs1oE[i]->Draw("same,p");
        DrawGammaSetMarkerTF1(fitResoF1oE[i], lineStylePart[i], 3, colorPart[i]-6);
        fitResoF1oE[i]->Draw("same");
        legendResPart1oE->AddEntry(histResoEVs1oE[i],Form("%s, r", labelPart[i].Data()),"p");
        legendResPart1oE->AddEntry(fitReso1oE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitReso1oE[i]->GetParameter(1),fitReso1oE[i]->GetParameter(0)),"l");
        legendResPart1oE->AddEntry(histResoFEVs1oE[i],Form("%s, l", labelPart[i].Data()),"p");
        legendResPart1oE->AddEntry(fitResoF1oE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoF1oE[i]->GetParameter(1),fitResoF1oE[i]->GetParameter(0)),"l");
      } else {      
        legendResPart1oE->AddEntry(histResoEVs1oE[i],Form("%s", labelPart[i].Data()),"p");
        legendResPart1oE->AddEntry(fitReso1oE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitReso1oE[i]->GetParameter(1),fitReso1oE[i]->GetParameter(0)),"l");
      }
    }
    legendResPart1oE->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso2->Print(Form("%s/Resolution_Fitted_1oE_AllPart.%s", outputDir.Data(), suffix.Data()));

  histResoDummy->Draw();
  
    TLegend* legendResPartHighest1oE  = GetAndSetLegend2(0.14, 0.95-(nActivePart*0.85*textSizeLabelsRel), 0.6, 0.95,0.85*textSizeLabelsPixel, 2, "", 43, 0.2);
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      DrawGammaSetMarker(histResoEVs1oEhighest[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histResoEVs1oEhighest[i]->Draw("same,p");
      DrawGammaSetMarkerTF1(fitReso1oEhighest[i], lineStylePart[i], 3, colorPart[i]+1);
      fitReso1oEhighest[i]->Draw("same");
      if (i < 2 && doFitting){
        DrawGammaSetMarker(histResoFEVs1oEhighest[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i]-6, colorPart[i]-6);
        histResoFEVs1oEhighest[i]->Draw("same,p");
        DrawGammaSetMarkerTF1(fitResoF1oEhighest[i], lineStylePart[i], 3, colorPart[i]-6);
        fitResoF1oEhighest[i]->Draw("same");
        legendResPartHighest1oE->AddEntry(histResoEVs1oEhighest[i],Form("%s, r", labelPart[i].Data()),"p");
        legendResPartHighest1oE->AddEntry(fitReso1oEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitReso1oEhighest[i]->GetParameter(1),fitReso1oEhighest[i]->GetParameter(0)),"l");
        legendResPartHighest1oE->AddEntry(histResoFEVs1oEhighest[i],Form("%s, l", labelPart[i].Data()),"p");
        legendResPartHighest1oE->AddEntry(fitResoF1oEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoF1oEhighest[i]->GetParameter(1),fitResoF1oEhighest[i]->GetParameter(0)),"l");
      } else {      
        legendResPartHighest1oE->AddEntry(histResoEVs1oEhighest[i],Form("%s", labelPart[i].Data()),"p");
        legendResPartHighest1oE->AddEntry(fitReso1oEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitReso1oEhighest[i]->GetParameter(1),fitReso1oEhighest[i]->GetParameter(0)),"l");
      }
    }
    legendResPartHighest1oE->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso2->Print(Form("%s/Resolution_FittedHighest_1oE_AllPart.%s", outputDir.Data(), suffix.Data()));
   
  histResoDummy->Draw();
    if (isEMCal){
      graphReq1oEEM->Draw("same,e3");
    } else { 
      graphReq1oEEM->Draw("same,e3");
      graphReq1oEHad->Draw("same,e3");
    }

    histResoDummy->Draw("axis,same");
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      histResoEVs1oEhighest[i]->Draw("same,p");
      fitReso1oEhighest[i]->Draw("same");
      if (i < 2 && doFitting){
        histResoFEVs1oEhighest[i]->Draw("same,p");
        fitResoF1oEhighest[i]->Draw("same");
      }
    }
    legendResPartHighest1oE->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (isEMCal)
      legendResYR1oEECal->Draw();
    else 
      legendResYR1oECAL->Draw();

  cReso2->Print(Form("%s/Resolution_FittedHighest_1oE_AllPart_wReq.%s", outputDir.Data(), suffix.Data()));

  for (Int_t etbin = minEtaBinCaloDis[dirCal]; etbin < maxEtaBinCaloDis[dirCal]; etbin++){
    histResoDummy->Draw();
    bool draw = 0;
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!useEta[i][etbin]) continue;
      if (i > 1 && isEMCal) continue;
      draw = 1;
      DrawGammaSetMarker(histResoEVs1oEEta[i][etbin], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histResoEVs1oEEta[i][etbin]->Draw("same,p");
      DrawGammaSetMarkerTF1(fitReso1oE[i], lineStylePart[i], 3, colorPart[i]+1);
      fitReso1oE[i]->Draw("same");
    }
    legendResPart->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%1.1f< #eta< %1.1f",partEtaCalo[etbin],partEtaCalo[etbin+1]),0.95,0.83,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
    
    if (draw)cReso2->Print(Form("%s/Resolution_Fitted_1oE_AllPart_Eta%d.%s", outputDir.Data(), etbin,suffix.Data()));
  }

  cReso->cd();
  DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
  TH2F * histResoDummyMean    = new TH2F("histResoDummyMean","histResoDummyMean",1000,0.3, 170,1000,0.2, 1.5);
  SetStyleHistoTH2ForGraphs(histResoDummyMean, "#it{E} (GeV)","#mu", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummyMean->GetXaxis()->SetNoExponent();
  histResoDummyMean->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummyMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummyMean->Draw();

    TLegend* legendResPartMu  = GetAndSetLegend2(0.13, 0.96-(nActivePart*0.85*textSizeLabelsRel), 0.43, 0.96,0.85*textSizeLabelsPixel, 2, "", 43, 0.22);
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      DrawGammaSetMarker(histMeanEVsE[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histMeanEVsE[i]->Draw("same,p");
      fitMeanE[i] = new TF1(Form("fit_mean_%sE", readNames[i].Data()), "[0]", meanFitEM[dirCal][0], meanFitEM[dirCal][1]);
      histMeanEVsE[i]->Fit(fitMeanE[i],"RMNE");
      fitMeanE[i]->SetNpx(10000);
      fitMeanE[i]->Draw("same");
      DrawGammaSetMarkerTF1(fitMeanE[i], 3, 3, colorPart[i]+1);
      legendResPartMu->AddEntry(histMeanEVsE[i],Form("%s", labelPart[i].Data()),"p");
      legendResPartMu->AddEntry(fitMeanE[i],Form("#mu =  %1.2f",fitMeanE[i]->GetParameter(0)),"l");
    }
    legendResPartMu->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/Mean_allPart.%s", outputDir.Data(), suffix.Data()));

  TFile* outputFile = new TFile(Form("%s/ResolutionPostProcessing.root", outputDir.Data()),"RECREATE");
  outputFile->cd();
  
  for (Int_t i = 0; i <nParticles-1; i++){
    if (!usePart[i]) continue;
    if (histSigmaEVsE[i]) histSigmaEVsE[i]->Write();
    if (histFWHMEVsE[i]) histFWHMEVsE[i]->Write();
    if (histResoEVsE[i]) histResoEVsE[i]->Write();
    if (fitResoE[i]) fitResoE[i]->Write();
    if (histResoFEVsE[i]) histResoFEVsE[i]->Write();
    if (fitResoFE[i]) fitResoFE[i]->Write();
    if (histMeanEVsE[i]) histMeanEVsE[i]->Write();
    if (fitMeanE[i]) fitMeanE[i]->Write();
    if (histResoEVs1oE[i]) histResoEVs1oE[i]->Write();
    if (fitReso1oE[i]) fitReso1oE[i]->Write();
    if (histResoFEVs1oE[i]) histResoFEVs1oE[i]->Write();
    if (fitResoF1oE[i]) fitResoF1oE[i]->Write();
    
    if (histSigmaEVsEhighest[i]) histSigmaEVsEhighest[i]->Write();
    if (histFWHMEVsEhighest[i]) histFWHMEVsEhighest[i]->Write();
    if (histResoEVsEhighest[i]) histResoEVsEhighest[i]->Write();
    if (fitResoEhighest[i]) fitResoEhighest[i]->Write();
    if (histResoFEVsEhighest[i]) histResoFEVsEhighest[i]->Write();
    if (fitResoFEhighest[i]) fitResoFEhighest[i]->Write();
    if (histMeanEVsEhighest[i]) histMeanEVsEhighest[i]->Write();
    if (histResoEVs1oEhighest[i]) histResoEVs1oEhighest[i]->Write();
    if (fitReso1oEhighest[i]) fitReso1oEhighest[i]->Write();
    if (histResoFEVs1oEhighest[i]) histResoFEVs1oEhighest[i]->Write();
    if (fitResoF1oEhighest[i]) fitResoF1oEhighest[i]->Write();
    
    for (Int_t ebin = 0; ebin < nEne; ebin++){
      if (!useE[i][ebin]) continue;
      if (hist1DResoEbins[i][ebin]) hist1DResoEbins[i][ebin]->Write();
      if (fitResoEbins[i][ebin]) fitResoEbins[i][ebin]->Write();
      if (hist1DResoEbinshighests[i][ebin]) hist1DResoEbinshighests[i][ebin]->Write();
      if (fitResoEbinshighest[i][ebin]) fitResoEbinshighest[i][ebin]->Write();
    }
    
    for (Int_t etbin = minEtaBinCaloDis[dirCal]; etbin < maxEtaBinCaloDis[dirCal]; etbin++){
      if (histSigmaEVsEEta[i][etbin]) histSigmaEVsEEta[i][etbin]->Write();
      if (histFWHMEVsEEta[i][etbin]) histFWHMEVsEEta[i][etbin]->Write();
      if (histResoEVsEEta[i][etbin]) histResoEVsEEta[i][etbin]->Write();
      if (fitResoEEta[i][etbin]) fitResoEEta[i][etbin]->Write();
      if (histResoFEVsEEta[i][etbin]) histResoFEVsEEta[i][etbin]->Write();
      if (fitResoFEEta[i][etbin]) fitResoFEEta[i][etbin]->Write();
      if (histMeanEVsEEta[i][etbin]) histMeanEVsEEta[i][etbin]->Write();
      if (histResoEVs1oEEta[i][etbin]) histResoEVs1oEEta[i][etbin]->Write();
      if (fitReso1oEEta[i][etbin]) fitReso1oEEta[i][etbin]->Write();
      if (histResoFEVs1oEEta[i][etbin]) histResoFEVs1oEEta[i][etbin]->Write();
      if (fitResoF1oEEta[i][etbin]) fitResoF1oEEta[i][etbin]->Write();
      for (Int_t ebin = 0; ebin < nEne; ebin++){
        if (!useEEta[i][ebin][etbin]) continue;
        if (hist1DResoEbinsEta[i][ebin][etbin]) hist1DResoEbinsEta[i][ebin][etbin]->Write();
        if (fitResoEbinsEta[i][ebin][etbin]) fitResoEbinsEta[i][ebin][etbin]->Write();
      }
    }
  }
  
  outputFile->Write();
  outputFile->Close();
  
}
