#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

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
  
  
  const static Double_t partE[]   = {   0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 
                                        5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 
                                        12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 50, 
                                        67.5, 75, 100, 150
                                    };

  TGraphErrors* graphReq          = new TGraphErrors(nEne); 
  TGraphErrors* graphReq1oE       = new TGraphErrors(nEne); 
  for (Int_t i= 0; i < nEne; i++){
//     cout << "ebin:" << partE[i] << endl;
    Double_t value[4]             = {0};
    value[0]  = sqrtETerm[0]/TMath::Sqrt(partE[i])+linTerm[0];
    value[1]  = sqrtETerm[1]/TMath::Sqrt(partE[i])+linTerm[0];
    value[2]  = sqrtETerm[0]/TMath::Sqrt(partE[i])+linTerm[1];
    value[3]  = sqrtETerm[1]/TMath::Sqrt(partE[i])+linTerm[1];
//     cout << value[0] << "\t" << value[1] << "\t"<< value[2] << "\t"<< value[3] << "\t"<< endl;
    Double_t max = -10000;
    Double_t min = 10000;
    for (Int_t k = 0; k < 4; k++){ 
        if (min > value[k]) min = value[k];
        if (max < value[k]) max = value[k];
    }
    graphReq->SetPoint(i, partE[i], (min+max)/2 );
    graphReq->SetPointError(i, 0.05, TMath::Abs(min-max)/2 );
    graphReq1oE->SetPoint(i, 1/TMath::Sqrt(partE[i]), (min+max)/2 );
    graphReq1oE->SetPointError(i, 0.005, TMath::Abs(min-max)/2 );
  }
  graphReq1oE->Sort();
  DrawGammaSetMarkerTGraphErr(  graphReq, 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
  DrawGammaSetMarkerTGraphErr(  graphReq1oE, 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
  TLegend* legendResYR    = GetAndSetLegend2(0.14, 0.14, 0.6, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYR->AddEntry(graphReq, "YR requirement","f");
  legendResYR->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETerm[0], sqrtETerm[1], linTerm[0], linTerm[1]),"");
  TLegend* legendResYR1oE = GetAndSetLegend2(0.58, 0.14, 0.97, 0.14+(2*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendResYR1oE->AddEntry(graphReq1oE, "YR requirement","f");
  legendResYR1oE->AddEntry((TObject*)0, Form("#sigma/#it{E} =  (%1.0f-%1.0f)/#sqrt{#it{E}} #oplus (%1.0f-%1.0f)", sqrtETerm[0], sqrtETerm[1], linTerm[0], linTerm[1]),"");
  
  
  Int_t useE[7][nEne]             = {{ 0}};
  Int_t usePart[7]                = { 0};
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
  Int_t nActiveEpart[7]             = {0};
  Int_t nActiveEtapart[7]           = {0};
  
  TFile* inputFile                  = new TFile(inputFileName.Data());
  TH2D* hist2DResoE[7]              = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH2D* hist2DResoEhighest[7]       = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1D* hist1DResoEbins[7][nEne]    = {{NULL}};
  TH1D* hist1DResoEbinshighests[7][nEne]    = {{NULL}};
  TH1F* histSigmaEVsE[7]            = {NULL};
  TH1F* histResoEVsE[7]             = {NULL};
  TH1F* histMeanEVsE[7]             = {NULL};
  TH1F* histResoEVs1oE[7]           = {NULL};
  TF1* fitMeanE[7]                  = {NULL};
  TF1* fitResoE[7]                  = {NULL};
  TF1* fitReso1oE[7]                = {NULL};
  TF1* fitResoEbins[7][nEne]        = {{NULL}};
  TH1F* histSigmaEVsEhighest[7]            = {NULL};
  TH1F* histResoEVsEhighest[7]             = {NULL};
  TH1F* histMeanEVsEhighest[7]             = {NULL};
  TH1F* histResoEVs1oEhighest[7]           = {NULL};
  TF1* fitResoEhighest[7]                  = {NULL};
  TF1* fitReso1oEhighest[7]                = {NULL};
  TF1* fitResoEbinshighest[7][nEne]        = {{NULL}};

  TH2D* hist2DResoEEta[7][nEta]              = {{NULL}};
  TH1D* hist1DResoEbinsEta[7][nEne][nEta]    = {{{NULL}}};
  TH1F* histSigmaEVsEEta[7][nEta]            = {{NULL}};
  TH1F* histResoEVsEEta[7][nEta]             = {{NULL}};
  TH1F* histMeanEVsEEta[7][nEta]             = {{NULL}};
  TH1F* histResoEVs1oEEta[7][nEta]           = {{NULL}};
  TF1* fitResoEEta[7][nEta]                  = {{NULL}};
  TF1* fitReso1oEEta[7][nEta]                = {{NULL}};
  TF1* fitResoEbinsEta[7][nEne][nEta]        = {{NULL}};
  Int_t useEEta[7][nEne][nEta]               = {{{ 0}}};
  Int_t useEta[7][nEta]                      = {{0}};
  
  for (Int_t i = 0; i < nParticles; i++){
    histSigmaEVsE[i]      = new TH1F(Form("h_%ssigma_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histMeanEVsE[i]       = new TH1F(Form("h_%smean_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVsE[i]       = new TH1F(Form("h_%sreso_E", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVs1oE[i]     = new TH1F(Form("h_%sreso_1oE", readNames[i].Data()), "", 1000, 0., 4.);
    histSigmaEVsEhighest[i]      = new TH1F(Form("h_%ssigma_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histMeanEVsEhighest[i]       = new TH1F(Form("h_%smean_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVsEhighest[i]       = new TH1F(Form("h_%sreso_Ehighest", readNames[i].Data()), "", 401, -0.25, 200.25);
    histResoEVs1oEhighest[i]     = new TH1F(Form("h_%sreso_1oEhighest", readNames[i].Data()), "", 1000, 0., 4.);
    
    
    TString tempName = Form("%s/h_CRH_EReso_%sE_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempName = Form("%s/h_CRH_EResoComb_%sE_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
//     cout << tempName << endl;
    TString tempName3 = Form("%s/h_CRH_EReso_%sEhighest_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    if (isComb) tempName3 = tempName;
//     cout << tempName << endl;
    hist2DResoE[i]    = (TH2D*)inputFile->Get(tempName.Data());
    hist2DResoEhighest[i]    = (TH2D*)inputFile->Get(tempName3.Data());
    if (hist2DResoE[i]->GetEntries() > 0){
      hist2DResoE[i]->Sumw2();
      hist2DResoEhighest[i]->Sumw2();
      for (int ebin = 0; ebin < nEne; ebin++){
        hist1DResoEbins[i][ebin]    = (TH1D*)hist2DResoE[i]->ProjectionY(Form("projection_dE_%s%d", readNames[i].Data(), ebin), 
                                                                        hist2DResoE[i]->GetXaxis()->FindBin(partE[ebin]-0.05), hist2DResoE[i]->GetXaxis()->FindBin(partE[ebin]+0.05),"e"); 
        if (hist1DResoEbins[i][ebin]->GetEntries() > 10 ){
          useE[i][ebin]         = 1;
          usePart[i]            = 1;
          nActiveEpart[i]++;
          
          Double_t mean     = FindLargestBin1DHist(hist1DResoEbins[i][ebin], 0.4, 2);
//           Double_t mean     = hist1DResoEbins[i][ebin]->GetMean();
          Double_t meanErr  = hist1DResoEbins[i][ebin]->GetMeanError();
//           Double_t sigma    = hist1DResoEbins[i][ebin]->GetStdDev();
//           Double_t sigmaErr = hist1DResoEbins[i][ebin]->GetStdDevError();
          Double_t sigma    = hist1DResoEbins[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoEbins[i][ebin]->GetRMSError();
          
          if (doFitting){
//             cout << labelPart[i] <<"\t" << partE[ebin] << "\t" << mean << "\t" << sigma<< endl;;
//             fitResoEbins[i][ebin]     = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))", mean-3*sigma, mean+2*sigma);
            fitResoEbins[i][ebin]     = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "gaus", mean-3*sigma, mean+2*sigma);
            fitResoEbins[i][ebin]->SetParLimits(1, 0.6*mean, 1.2*mean);
            fitResoEbins[i][ebin]->SetParameter(1, mean);
            fitResoEbins[i][ebin]->SetParLimits(2, 0.6*sigma,sigma);
            fitResoEbins[i][ebin]->SetParameter(2, sigma);
//             fitResoEbins[i][ebin]->SetParLimits(0, 0.9*hist1DResoEbins[i][ebin]->GetMaximum(), 10*hist1DResoEbins[i][ebin]->GetMaximum());
//             fitResoEbins[i][ebin]->SetParameter(0, hist1DResoEbins[i][ebin]->GetMaximum());
            hist1DResoEbins[i][ebin]->Rebin(rebinRes);
            hist1DResoEbins[i][ebin]->Fit(fitResoEbins[i][ebin],"L0RMEQ","",mean-1*sigma, mean+1.2*sigma);
//             hist1DResoEbins[i][ebin]->Fit(fitResoEbins[i][ebin],"L0RME", "", fitResoEbins[i][ebin]->GetParameter(1)-0.5*fitResoEbins[i][ebin]->GetParameter(2), fitResoEbins[i][ebin]->GetParameter(1)+1*fitResoEbins[i][ebin]->GetParameter(2));
            
            mean     = fitResoEbins[i][ebin]->GetParameter(1);
            meanErr  = fitResoEbins[i][ebin]->GetParError(1);
            sigma    = fitResoEbins[i][ebin]->GetParameter(2);
            sigmaErr = fitResoEbins[i][ebin]->GetParError(2);
          }
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));

          
          histResoEVsE[i]->SetBinContent(histResoEVsE[i]->FindBin(partE[ebin]),  100*reso);
          histResoEVsE[i]->SetBinError(histResoEVsE[i]->FindBin(partE[ebin]), resoErr);
          histResoEVs1oE[i]->SetBinContent(histResoEVs1oE[i]->FindBin(1./TMath::Sqrt(partE[ebin])),100*reso);
          histResoEVs1oE[i]->SetBinError(histResoEVs1oE[i]->FindBin(1./TMath::Sqrt(partE[ebin])),resoErr);

          histMeanEVsE[i]->SetBinContent(histResoEVsE[i]->FindBin(partE[ebin]), mean);
          histMeanEVsE[i]->SetBinError(histResoEVsE[i]->FindBin(partE[ebin]), meanErr);
          histSigmaEVsE[i]->SetBinContent(histResoEVsE[i]->FindBin(partE[ebin]), sigma);
          histSigmaEVsE[i]->SetBinError(histResoEVsE[i]->FindBin(partE[ebin]), sigmaErr);
        } else {
          useE[i][ebin] = 0;
        }
      }
      for (int ebin = 0; ebin < nEne; ebin++){
        hist1DResoEbinshighests[i][ebin]    = (TH1D*)hist2DResoEhighest[i]->ProjectionY(Form("projection_dEhighest_%s%d", readNames[i].Data(), ebin), 
                                                                        hist2DResoEhighest[i]->GetXaxis()->FindBin(partE[ebin]-0.05), hist2DResoEhighest[i]->GetXaxis()->FindBin(partE[ebin]+0.05),"e"); 
        if (hist1DResoEbinshighests[i][ebin]->GetEntries() > 10 ){
          
          Double_t mean     = FindLargestBin1DHist(hist1DResoEbinshighests[i][ebin], 0.4, 2);
//           Double_t mean     = hist1DResoEbins[i][ebin]->GetMean();
          Double_t meanErr  = hist1DResoEbinshighests[i][ebin]->GetMeanError();
//           Double_t sigma    = hist1DResoEbins[i][ebin]->GetStdDev();
//           Double_t sigmaErr = hist1DResoEbins[i][ebin]->GetStdDevError();
          Double_t sigma    = hist1DResoEbinshighests[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoEbinshighests[i][ebin]->GetRMSError();
          
          if (doFitting){
//             cout << labelPart[i] <<"\t" << partE[ebin] << "\t" << mean << "\t" << sigma<< endl;;
//             fitResoEbins[i][ebin]     = new TF1(Form("fit_reso_%s%d", readNames[i].Data(), ebin), "(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))", mean-3*sigma, mean+2*sigma);
            fitResoEbinshighest[i][ebin]     = new TF1(Form("fit_resohighest_%s%d", readNames[i].Data(), ebin), "gaus", mean-3*sigma, mean+2*sigma);
            fitResoEbinshighest[i][ebin]->SetParLimits(1, 0.6*mean, 1.2*mean);
            fitResoEbinshighest[i][ebin]->SetParameter(1, mean);
            fitResoEbinshighest[i][ebin]->SetParLimits(2, 0.6*sigma,sigma);
            fitResoEbinshighest[i][ebin]->SetParameter(2, sigma);
//             fitResoEbins[i][ebin]->SetParLimits(0, 0.9*hist1DResoEbins[i][ebin]->GetMaximum(), 10*hist1DResoEbins[i][ebin]->GetMaximum());
//             fitResoEbins[i][ebin]->SetParameter(0, hist1DResoEbins[i][ebin]->GetMaximum());
            hist1DResoEbinshighests[i][ebin]->Rebin(rebinRes);
            hist1DResoEbinshighests[i][ebin]->Fit(fitResoEbinshighest[i][ebin],"L0RMEQ","",mean-1*sigma, mean+1.2*sigma);
//             hist1DResoEbins[i][ebin]->Fit(fitResoEbins[i][ebin],"L0RME", "", fitResoEbins[i][ebin]->GetParameter(1)-0.5*fitResoEbins[i][ebin]->GetParameter(2), fitResoEbins[i][ebin]->GetParameter(1)+1*fitResoEbins[i][ebin]->GetParameter(2));
            
            mean     = fitResoEbinshighest[i][ebin]->GetParameter(1);
            meanErr  = fitResoEbinshighest[i][ebin]->GetParError(1);
            sigma    = fitResoEbinshighest[i][ebin]->GetParameter(2);
            sigmaErr = fitResoEbinshighest[i][ebin]->GetParError(2);
          }
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));

          
          histResoEVsEhighest[i]->SetBinContent(histResoEVsEhighest[i]->FindBin(partE[ebin]),  100*reso);
          histResoEVsEhighest[i]->SetBinError(histResoEVsEhighest[i]->FindBin(partE[ebin]), resoErr);
          histResoEVs1oEhighest[i]->SetBinContent(histResoEVs1oEhighest[i]->FindBin(1./TMath::Sqrt(partE[ebin])),100*reso);
          histResoEVs1oEhighest[i]->SetBinError(histResoEVs1oEhighest[i]->FindBin(1./TMath::Sqrt(partE[ebin])),resoErr);

          histMeanEVsEhighest[i]->SetBinContent(histResoEVsEhighest[i]->FindBin(partE[ebin]), mean);
          histMeanEVsEhighest[i]->SetBinError(histResoEVsEhighest[i]->FindBin(partE[ebin]), meanErr);
          histSigmaEVsEhighest[i]->SetBinContent(histResoEVsEhighest[i]->FindBin(partE[ebin]), sigma);
          histSigmaEVsEhighest[i]->SetBinError(histResoEVsEhighest[i]->FindBin(partE[ebin]), sigmaErr);
        } else {
          useE[i][ebin] = 0;
        }
      }
      for (Int_t etbin = 0; etbin < nEta; etbin++){
        histSigmaEVsEEta[i][etbin]      = new TH1F(Form("h_%ssigma_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histMeanEVsEEta[i][etbin]       = new TH1F(Form("h_%smean_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histResoEVsEEta[i][etbin]       = new TH1F(Form("h_%sreso_E_Eta_%d", readNames[i].Data(), etbin), "", 401, -0.25, 200.25);
        histResoEVs1oEEta[i][etbin]     = new TH1F(Form("h_%sreso_1oE_Eta_%d", readNames[i].Data(), etbin), "", 1000, 0., 4.);
        TString tempName2 = Form("%s/h_CRH_EReso_%sE_Eta_%s_%s_%d",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data(), etbin );
        if (isComb) tempName2 = Form("%s/h_CRH_EResoComb_%sE_Eta_%s_%s_%d",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data(), etbin );
//         cout << tempName2 << endl;

        hist2DResoEEta[i][etbin]    = (TH2D*)inputFile->Get(tempName2.Data());
        if (!hist2DResoEEta[i][etbin]) continue;
        if (hist2DResoEEta[i][etbin]->GetEntries() > 0){
          hist2DResoEEta[i][etbin]->Sumw2();
          for (int ebin = 0; ebin < nEne; ebin++){
            hist1DResoEbinsEta[i][ebin][etbin]    = (TH1D*)hist2DResoEEta[i][etbin]->ProjectionY(Form("projection_dE_%s%d_%d", readNames[i].Data(), ebin, etbin), 
                                                                            hist2DResoEEta[i][etbin]->GetXaxis()->FindBin(partE[ebin]-0.05), hist2DResoEEta[i][etbin]->GetXaxis()->FindBin(partE[ebin]+0.05),"e"); 
            if (hist1DResoEbinsEta[i][ebin][etbin]->GetEntries() > 10){
              useEEta[i][ebin][etbin]  = 1;
              if (!useEta[i][etbin]) nActiveEtapart[i]++;
              useEta[i][etbin]         = 1;
              
              Double_t meanEta     = FindLargestBin1DHist(hist1DResoEbinsEta[i][ebin][etbin], 0.4, 2);
//               Double_t meanEta     = hist1DResoEbinsEta[i][ebin][etbin]->GetMean();
              Double_t meanEtaErr  = hist1DResoEbinsEta[i][ebin][etbin]->GetMeanError();
//               Double_t sigmaEta    = hist1DResoEbinsEta[i][ebin][etbin]->GetStdDev();
//               Double_t sigmaEtaErr = hist1DResoEbinsEta[i][ebin][etbin]->GetStdDevError();
              Double_t sigmaEta    = hist1DResoEbinsEta[i][ebin][etbin]->GetRMS();
              Double_t sigmaEtaErr = hist1DResoEbinsEta[i][ebin][etbin]->GetRMSError();

              if (doFitting){
                fitResoEbinsEta[i][ebin][etbin]     = new TF1(Form("fit_reso_%s%d_%d", readNames[i].Data(), ebin, etbin), "gaus", meanEta-3*sigmaEta, meanEta+2*sigmaEta);
//                 fitResoEbinsEta[i][ebin][etbin]     = new TF1(Form("fit_reso_%s%d_%d", readNames[i].Data(), ebin, etbin), "(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))", meanEta-3*sigmaEta, meanEta+2*sigmaEta);
                 
                fitResoEbinsEta[i][ebin][etbin]->SetParLimits(1, 0.6*meanEta, 1.2*meanEta);
                fitResoEbinsEta[i][ebin][etbin]->SetParameter(1, meanEta);
                fitResoEbinsEta[i][ebin][etbin]->SetParLimits(2, 0.6*sigmaEta, sigmaEta);
                fitResoEbinsEta[i][ebin][etbin]->SetParameter(2, sigmaEta);
//                 fitResoEbinsEta[i][ebin][etbin]->SetParLimits(0, 0.9*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 10*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum());
//                 fitResoEbinsEta[i][ebin][etbin]->SetParameter(0, hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum());
                fitResoEbinsEta[i][ebin][etbin]->SetParLimits(1, 0.6*fitResoEbins[i][ebin]->GetParameter(1), 1.4*fitResoEbins[i][ebin]->GetParameter(1));
                fitResoEbinsEta[i][ebin][etbin]->SetParameter(1, fitResoEbins[i][ebin]->GetParameter(1));
                fitResoEbinsEta[i][ebin][etbin]->SetParLimits(2, 0.6*fitResoEbins[i][ebin]->GetParameter(2), 1.4*fitResoEbins[i][ebin]->GetParameter(2));
                fitResoEbinsEta[i][ebin][etbin]->SetParameter(2, fitResoEbins[i][ebin]->GetParameter(2));                
                hist1DResoEbinsEta[i][ebin][etbin]->Rebin(rebinRes);
                hist1DResoEbinsEta[i][ebin][etbin]->Fit(fitResoEbinsEta[i][ebin][etbin],"L0RMEQ","", meanEta-1*sigmaEta, meanEta+1.2*sigmaEta);
//                 hist1DResoEbinsEta[i][ebin][etbin]->Fit(fitResoEbinsEta[i][ebin][etbin],"L0RMEQ", "", fitResoEbinsEta[i][ebin][etbin]->GetParameter(1)-1*fitResoEbinsEta[i][ebin][etbin]->GetParameter(2), fitResoEbinsEta[i][ebin][etbin]->GetParameter(1)+1.5*fitResoEbinsEta[i][ebin][etbin]->GetParameter(2));

                meanEta     = fitResoEbinsEta[i][ebin][etbin]->GetParameter(1);
                meanEtaErr  = fitResoEbinsEta[i][ebin][etbin]->GetParError(1);
                sigmaEta    = fitResoEbinsEta[i][ebin][etbin]->GetParameter(2);
                sigmaEtaErr = fitResoEbinsEta[i][ebin][etbin]->GetParError(2);
              }
              Double_t resoEta     = sigmaEta/meanEta;
              Double_t resoEtaErr  = TMath::Sqrt(TMath::Power(sigmaEtaErr/sigmaEta,2)+TMath::Power(meanEtaErr/meanEta,2));

              
              histResoEVsEEta[i][etbin]->SetBinContent(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), 100*resoEta);
              histResoEVsEEta[i][etbin]->SetBinError(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), resoEtaErr);
              histResoEVs1oEEta[i][etbin]->SetBinContent(histResoEVs1oEEta[i][etbin]->FindBin(1./TMath::Sqrt(partE[ebin])), 100*resoEta);
              histResoEVs1oEEta[i][etbin]->SetBinError(histResoEVs1oEEta[i][etbin]->FindBin(1./TMath::Sqrt(partE[ebin])), resoEtaErr);

              histMeanEVsEEta[i][etbin]->SetBinContent(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), meanEta);
              histMeanEVsEEta[i][etbin]->SetBinError(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), meanEtaErr);
              histSigmaEVsEEta[i][etbin]->SetBinContent(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), sigmaEta);
              histSigmaEVsEEta[i][etbin]->SetBinError(histResoEVsEEta[i][etbin]->FindBin(partE[ebin]), sigmaEtaErr);
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
  
   // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  for (Int_t i = 0; i < nParticles-1; i++){
    if (!usePart[i]) continue;
    TCanvas* cPNG;
    split_canvas(cPNG, Form("cPNG%d",i), nActiveEpart[i]);
    Int_t padnum =0;
    std::cout << labelPart[i] << std::endl;
      
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
      DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin)-histSigmaEVsE[i]->GetBinContent(bin), 
                     histMeanEVsE[i]->GetBinContent(bin)-histSigmaEVsE[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin)+histSigmaEVsE[i]->GetBinContent(bin), 
                     histMeanEVsE[i]->GetBinContent(bin)+histSigmaEVsE[i]->GetBinContent(bin), 
                     0., 0.05*hist1DResoEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNG->Print(Form("%s/Energy_%sResolution.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));

    for (Int_t etbin = 0; etbin < nEta; etbin++){
      if (!useEta[i][etbin]) continue;
      Double_t etaMinCurr = partEta[etbin];
      Double_t etaMaxCurr = partEta[etbin+1];
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
        Int_t bin = histMeanEVsE[i]->FindBin(partE[ebin]);
        DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin), 
                      histMeanEVsE[i]->GetBinContent(bin), 
                      0., 0.05*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 1, kGray+2, 2);
        DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin)-histSigmaEVsE[i]->GetBinContent(bin), 
                      histMeanEVsE[i]->GetBinContent(bin)-histSigmaEVsE[i]->GetBinContent(bin), 
                      0., 0.05*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 1, kRed+2, 2);
        DrawGammaLines(histMeanEVsE[i]->GetBinContent(bin)+histSigmaEVsE[i]->GetBinContent(bin), 
                      histMeanEVsE[i]->GetBinContent(bin)+histSigmaEVsE[i]->GetBinContent(bin), 
                      0., 0.05*hist1DResoEbinsEta[i][ebin][etbin]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNG->Print(Form("%s/Energy_%sResolution_Etabin%d.%s", outputDir.Data(), readNames[i].Data(), etbin, suffix.Data()));

    }
    
    if (maxFitRes1oE ==  -1) maxFitRes1oE  = 1/TMath::Sqrt(minEPart[i]);
    if (minFitRes1oE ==  -1) minFitRes1oE  = 1/TMath::Sqrt(maxEPart[i]);
    if (maxFitResE ==  -1) maxFitResE    = maxEPart[i];
    if (minFitResE ==  -1) minFitResE    = minEPart[i];
    if (maxFitMeanE ==  -1) maxFitMeanE  = maxEPart[i];
    if (minFitMeanE ==  -1) minFitMeanE  = minEPart[i];

    // RESOLUTION PLOT
    TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
    cReso->SetLogx(1);
    DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
    TH2F * histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0.3, 170,1000,0, maxResoPer);
    SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    histResoDummy->GetXaxis()->SetNoExponent();
    histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    histResoDummy->Draw();

    TLegend* legendResLabE  = GetAndSetLegend2(0.14, 0.95-0.85*textSizeLabelsRel, 0.35, 0.95,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);

    DrawGammaSetMarker(histResoEVsE[i], markerStylePart[i], markerSizePart[i], colorPart[i], colorPart[i]);
    histResoEVsE[i]->Draw("same,p");
//     fitResoE[i] = new TF1(Form("fit_reso_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x)))", minEPart[i],maxEPart[i]);
    fitResoE[i] = new TF1(Form("fit_reso_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x))+([2]/x)", minFitResE, maxFitResE);
    fitResoE[i]->SetParLimits(0,0.,20);
    fitResoE[i]->SetParLimits(1,0.,90);
    fitResoE[i]->SetParameter(2,0.01);
    fitResoE[i]->SetParLimits(2,0.,0.3);
    histResoEVsE[i]->Fit(fitResoE[i],"RMNE");
    fitResoE[i]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fitResoE[i], 3, 3, colorPart[i]+1);
    fitResoE[i]->SetRange(minEPart[i], maxEPart[i]);
    fitResoE[i]->Draw("same");
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendResLabE->AddEntry(fitResoE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
    legendResLabE->Draw();

    cReso->Print(Form("%s/Resolution_%sFitted.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));

    DrawGammaSetMarker(histResoEVsEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]+2, colorPart[i]+2);
    histResoEVsEhighest[i]->Draw("same,p");
//     fitResoE[i] = new TF1(Form("fit_reso_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x))+([2]/x)", minEPart[i],maxEPart[i]);
    fitResoEhighest[i] = new TF1(Form("fit_reso_%sE", readNames[i].Data()), "[0]+([1]/TMath::Sqrt(x))+([2]/x)", minFitResE,maxFitResE);
    fitResoEhighest[i]->SetParLimits(0,0.,20);
    fitResoEhighest[i]->SetParLimits(1,0.,90);
    fitResoEhighest[i]->SetParameter(2,0.01);
    fitResoEhighest[i]->SetParLimits(2,0.,0.3);
    histResoEVsEhighest[i]->Fit(fitResoEhighest[i],"RMNE");
    DrawGammaSetMarkerTF1(fitResoEhighest[i], 7, 3, colorPart[i]+2);
    fitResoEhighest[i]->SetRange(minEPart[i], maxEPart[i]);
    fitResoEhighest[i]->Draw("same");
    
    legendResLabE  = GetAndSetLegend2(0.14, 0.95-2*0.85*textSizeLabelsRel, 0.35, 0.95,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);
    legendResLabE->AddEntry(fitResoE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
    legendResLabE->AddEntry(fitResoEhighest[i], Form("%s highest: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitResoEhighest[i]->GetParameter(1),fitResoEhighest[i]->GetParameter(0)),"l");
    legendResLabE->Draw();

    
    cReso->Print(Form("%s/Resolution_%sFittedWithHighest.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
    histResoDummy->Draw();
    TLegend* legendResEEta  = GetAndSetLegend2(0.14, 0.95-((nActiveEtapart[i]+1)*0.85*textSizeLabelsRel), 0.35, 0.95,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);

    graphReq->Draw("same,e3");
    histResoDummy->Draw("axis,same");
    for (Int_t etbin = 0; etbin < nEta; etbin++){
      if (!useEta[i][etbin]) continue;
      DrawGammaSetMarker(histResoEVsEEta[i][etbin], markerStyleEta[etbin], markerSizeEta[etbin]*1.5, colorEta[etbin], colorEta[etbin]);
      histResoEVsEEta[i][etbin]->Draw("same,p");
      legendResEEta->AddEntry(histResoEVsEEta[i][etbin],Form("%1.1f< #eta< %1.1f", partEta[etbin], partEta[etbin+1]),"p");
    }
    DrawGammaSetMarkerTF1(fitResoE[i], 3, 3, kBlack);
    fitResoE[i]->Draw("same");
    legendResEEta->AddEntry(fitResoE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
    legendResYR->Draw();
    
    legendResEEta->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/Resolution_%sFitted_etaDep.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
      
    
  // RESOLUTION PLOT 2

    TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
    DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
    histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 1/TMath::Sqrt(minEPart[i])+0.1,1000,0, maxResoPer);
    SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    histResoDummy->GetXaxis()->SetNoExponent();
    histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    histResoDummy->Draw();
    
    TLegend* legendResLab1oE  = GetAndSetLegend2(0.14, 0.91-0.85*textSizeLabelsRel, 0.35, 0.91,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);

    DrawGammaSetMarker(histResoEVs1oE[i], markerStylePart[i], markerSizePart[i], colorPart[i], colorPart[i]);
    histResoEVs1oE[i]->Draw("same,p");
    fitReso1oE[i] = new TF1(Form("fit_reso_%s1oE", readNames[i].Data()), "[0]+[1]*x+[2]*TMath::Power(x,2)", minFitRes1oE,maxFitRes1oE); //+[2]*TMath::Power(x,2)
    fitReso1oE[i]->SetParameter(2,0.01);
    fitReso1oE[i]->SetParLimits(2,0.001, 0.3);

    histResoEVs1oE[i]->Fit(fitReso1oE[i],"RMNE");
    fitReso1oE[i]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fitReso1oE[i], 3, 3, colorPart[i]+1);
    fitReso1oE[i]->SetRange(1/TMath::Sqrt(maxEPart[i]), 1/TMath::Sqrt(minEPart[i]));
    fitReso1oE[i]->Draw("same");
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.14,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    legendResLab1oE->AddEntry(fitReso1oE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitReso1oE[i]->GetParameter(1),fitReso1oE[i]->GetParameter(0)),"l");
    legendResLab1oE->Draw();

    cReso2->Print(Form("%s/Resolution_%sFitted_1oE.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));


    DrawGammaSetMarker(histResoEVs1oEhighest[i], markerStylePart[i], markerSizePart[i], colorPart[i]+2, colorPart[i]+2);
    histResoEVs1oEhighest[i]->Draw("same,p");
    fitReso1oEhighest[i] = new TF1(Form("fit_reso_%s1oE", readNames[i].Data()), "[0]+[1]*x+[2]*TMath::Power(x,2)", minFitRes1oE,maxFitRes1oE); //+[2]*TMath::Power(x,2)
    fitReso1oEhighest[i]->SetParameter(2,0.01);
    fitReso1oEhighest[i]->SetParLimits(2,0.0, 0.3);
    histResoEVs1oEhighest[i]->Fit(fitReso1oEhighest[i],"RMNE");
    DrawGammaSetMarkerTF1(fitReso1oEhighest[i], 7, 3, colorPart[i]+2);
    fitReso1oEhighest[i]->SetRange(1/TMath::Sqrt(maxEPart[i]), 1/TMath::Sqrt(minEPart[i]));
    fitReso1oEhighest[i]->Draw("same");
    
    legendResLab1oE  = GetAndSetLegend2(0.14, 0.91-2*0.85*textSizeLabelsRel, 0.35, 0.91,0.85*textSizeLabelsPixel, 1, "", 43, 0.2);
    legendResLab1oE->AddEntry(fitReso1oE[i], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitReso1oE[i]->GetParameter(1),fitReso1oE[i]->GetParameter(0)),"l");
    legendResLab1oE->AddEntry(fitReso1oEhighest[i], Form("%s highest: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[i].Data(), fitReso1oEhighest[i]->GetParameter(1),fitReso1oEhighest[i]->GetParameter(0)),"l");
    legendResLab1oE->Draw();

    
    cReso2->Print(Form("%s/Resolution_%sFittedWithHighest_1oE.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
    
    histResoDummy->Draw();
    graphReq1oE->Draw("same,e3");
    histResoDummy->Draw("axis,same");
    for (Int_t etbin = 0; etbin < nEta; etbin++){
      if (!useEta[i][etbin]) continue;
      DrawGammaSetMarker(histResoEVs1oEEta[i][etbin], markerStyleEta[etbin], markerSizeEta[etbin]*1.5, colorEta[etbin], colorEta[etbin]);
      histResoEVs1oEEta[i][etbin]->Draw("same,p");
    }
    DrawGammaSetMarkerTF1(fitReso1oE[i], 3, 3, kBlack);
    fitReso1oE[i]->Draw("same");
    legendResYR1oE->Draw();
        
    legendResEEta->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    cReso2->Print(Form("%s/Resolution_%sFitted_1oE_etaDep.%s", outputDir.Data(), readNames[i].Data(), suffix.Data()));
  }
  
  TCanvas* cExampleBin = new TCanvas("cExampleBin","",0,0,1100,1000);
  DrawGammaCanvasSettings( cExampleBin, 0.095, 0.01, 0.01, 0.105);
  
    Float_t maximYRange = 0;
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!useE[i][exEbin]) continue;
      if (maximYRange < hist1DResoEbinshighests[i][exEbin]->GetMaximum())
        maximYRange = hist1DResoEbinshighests[i][exEbin]->GetMaximum();
    }
  
  
    TLegend* legendExamBin  = GetAndSetLegend2(0.14, 0.81-((nActivePart+1)/2*0.85*textSizeLabelsRel), 0.35, 0.81,0.85*textSizeLabelsPixel, 2, "", 43, 0.35);
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
    drawLatexAdd(perfLabel,0.14,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.14,0.87,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("%1.1f< #eta< %1.1f, #it{E} <%1.1f GeV", minEtaMax, maxEtaMax, partE[exEbin]),0.14,0.82,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    
  cExampleBin->Print(Form("%s/ExampleBin.%s", outputDir.Data(), suffix.Data()));
  
  
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
  cReso->SetLogx(1);
  DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
  TH2F * histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0.3, 170,1000,0, maxResoPer);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();

    if (isEMCal && nActivePart >2 ) nActivePart = 2;
    else nActivePart = nActivePart-1;
    TLegend* legendResPart  = GetAndSetLegend2(0.58, 0.86-(nActivePart*0.85*textSizeLabelsRel), 0.97, 0.86,0.85*textSizeLabelsPixel, 2, "", 43, 0.22);
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      DrawGammaSetMarker(histResoEVsE[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histResoEVsE[i]->Draw("same,p");
      DrawGammaSetMarkerTF1(fitResoE[i], lineStylePart[i], 3, colorPart[i]+1);
      fitResoE[i]->Draw("same");
      legendResPart->AddEntry(histResoEVsE[i],Form("%s", labelPart[i].Data()),"p");
      legendResPart->AddEntry(fitResoE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoE[i]->GetParameter(1),fitResoE[i]->GetParameter(0)),"l");
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
      legendResPartHighest->AddEntry(histResoEVsE[i],Form("%s", labelPart[i].Data()),"p");
      legendResPartHighest->AddEntry(fitResoEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitResoEhighest[i]->GetParameter(1),fitResoEhighest[i]->GetParameter(0)),"l");
    }
    legendResPartHighest->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/Resolution_FittedHighest_allPart.%s", outputDir.Data(), suffix.Data()));

  histResoDummy->Draw();
    graphReq->Draw("same,e3");
    histResoDummy->Draw("axis,same");
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      histResoEVsEhighest[i]->Draw("same,p");
      fitResoEhighest[i]->Draw("same");
    }
    legendResPartHighest->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendResYR->Draw();
    
  cReso->Print(Form("%s/Resolution_FittedHighest_allPart_wReq.%s", outputDir.Data(), suffix.Data()));

  
// RESOLUTION PLOT 2

  TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 1.5,1000,0, maxResoPer);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
  
    TLegend* legendResPart1oE  = GetAndSetLegend2(0.14, 0.95-(nActivePart*0.85*textSizeLabelsRel), 0.6, 0.95,0.85*textSizeLabelsPixel, 2, "", 43, 0.2);
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      DrawGammaSetMarker(histResoEVs1oE[i], markerStylePart[i], markerSizePart[i]*0.8, colorPart[i], colorPart[i]);
      histResoEVs1oE[i]->Draw("same,p");
      DrawGammaSetMarkerTF1(fitReso1oE[i], lineStylePart[i], 3, colorPart[i]+1);
      fitReso1oE[i]->Draw("same");
      legendResPart1oE->AddEntry(histResoEVs1oE[i],Form("%s", labelPart[i].Data()),"p");
      legendResPart1oE->AddEntry(fitReso1oE[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitReso1oE[i]->GetParameter(1),fitReso1oE[i]->GetParameter(0)),"l");
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
      legendResPartHighest1oE->AddEntry(histResoEVs1oEhighest[i],Form("%s", labelPart[i].Data()),"p");
      legendResPartHighest1oE->AddEntry(fitReso1oEhighest[i],Form("#sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",fitReso1oEhighest[i]->GetParameter(1),fitReso1oEhighest[i]->GetParameter(0)),"l");

    }
    legendResPartHighest1oE->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso2->Print(Form("%s/Resolution_FittedHighest_1oE_AllPart.%s", outputDir.Data(), suffix.Data()));
   
  histResoDummy->Draw();
    graphReq1oE->Draw("same,e3");
    histResoDummy->Draw("axis,same");
    for (Int_t i = 0; i <nParticles-1; i++){
      if (!usePart[i]) continue;
      if (i > 1 && isEMCal) continue;
      histResoEVs1oEhighest[i]->Draw("same,p");
      fitReso1oEhighest[i]->Draw("same");
    }
    legendResPartHighest1oE->Draw();
    drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendResYR1oE->Draw();
  cReso2->Print(Form("%s/Resolution_FittedHighest_1oE_AllPart_wReq.%s", outputDir.Data(), suffix.Data()));

  for (Int_t etbin = 0; etbin < nEta; etbin++){
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
    drawLatexAdd(Form("%1.1f< #eta< %1.1f",partEta[etbin],partEta[etbin+1]),0.95,0.83,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
    
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
      fitMeanE[i] = new TF1(Form("fit_reso_%sE", readNames[i].Data()), "[0]", minFitMeanE, maxFitMeanE);
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

  
}
