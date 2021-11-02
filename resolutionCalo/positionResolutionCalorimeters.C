#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void positionResolutionCalorimeters(
    TString inputFileName     = "output_ERH.root",
    TString caloName          = "BECAL",
    TString clusterizerName   = "MA",
    TString suffix            = "pdf",
    Bool_t doFitting          = kFALSE
){  

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  SetPlotStyle();
  TString collisionSystem           = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString perfLabel                 = "ECCE G4 simulation";
  
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;

  Double_t minEtaMax  = -4;
  Double_t maxEtaMax  = 4;
  Double_t minPhiMax  = -TMath::Pi()-0.1;
  Double_t maxPhiMax  = TMath::Pi()+0.1;
  Double_t maxReso    = 1.25;
  TString detLabel    = "";
  Bool_t isEMCal      = 0;
  Int_t rebinRes      = 1;
  TString caloNameRead = caloName;
  Bool_t isComb         = 0;
  
  Int_t exEbin          = 15;
  Int_t nbins           = 100; 
  Int_t region          = 0;
  
  if (caloName.CompareTo("EEMC") == 0){
    detLabel        = "EEMC (PbWO_{4} crystal)";
    minEtaMax       = -4;
    maxEtaMax       = -1.7;
    isEMCal         = 1;
    region          = 0;
  } else if (caloName.CompareTo("FEMC") == 0){
    detLabel        = "FEMC";
    minEtaMax       = 1.1;
    maxEtaMax       = 4;
    isEMCal         = 1;
    region          = 2;
  } else if (caloName.CompareTo("BECAL") == 0){
    detLabel        = "BECAL (Sci-glass)";
    minEtaMax       = -1.8;
    maxEtaMax       = 1.3;
    isEMCal         = 1;
    region          = 1;
  } else if (caloName.CompareTo("LFHCAL") == 0){
    detLabel        = "LFHCAL";
    minEtaMax       = 1.1;
    maxEtaMax       = 4;
    rebinRes        = 2;
    region          = 2;
  } else if (caloName.CompareTo("LFHCAL-wMat") == 0){
    detLabel        = "LFHCAL, w/ mat.";
    caloNameRead    = "LFHCAL";
    minEtaMax       = 1.1;
    maxEtaMax       = 4;
    rebinRes        = 2;
    region          = 2;
  } else if (caloName.CompareTo("LFHCAL-oFEMC") == 0){
    detLabel        = "LFHCAL, FEMC infront";
    caloNameRead    = "LFHCAL";
    minEtaMax       = 1.1;
    maxEtaMax       = 4;
    rebinRes        = 2;
    region          = 2;
  } else if (caloName.CompareTo("EHCAL") == 0){
    detLabel        = "EHCAL (STAR re-use)";
    minEtaMax       = -4;
    maxEtaMax       = -1.7;
    region          = 0;
  } else if (caloName.CompareTo("EHCAL-wMat") == 0){
    detLabel        = "EHCAL, w/ mat";
    caloNameRead    = "EHCAL";
    minEtaMax       = -4;
    maxEtaMax       = -1.7;
    region          = 0;
  } else if (caloName.CompareTo("CHCAL") == 0){
    detLabel        = "oHCAL, iHCal infront";
    caloNameRead    = "HCALOUT";
    minEtaMax       = -1.8;
    maxEtaMax       = 1.2;
    isComb          = 0;
    region          = 1;
  } else if (caloName.CompareTo("CHCAL-wMat") == 0){
    detLabel        = "oHCAL, w/ mat";
    caloNameRead    = "HCALOUT";
    minEtaMax       = -1.8;
    maxEtaMax       = 1.2;
    isComb          = 0;
    region          = 1;
  } else if (caloName.CompareTo("HCALIN") == 0){
    detLabel  = "iHCal";
    caloNameRead    = "HCALIN";
    minEtaMax       = -1.8;
    maxEtaMax       = 1.2;
    isComb          = 0;
    region          = 1;
  } else {
    minEtaMax       = -1.8;
    maxEtaMax       = 1.2;
  }
  // determine eta bins
  Int_t binEtaMin = 0;
  Int_t binEtaMax = 15;
  
  TString outputDir                 = Form("plotsResolution%s",caloName.Data());
  if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);


  //************************** Read data **************************************************
  SetPhiBins();
  
  TGraphErrors* graphReq          = new TGraphErrors(nP); 
  TGraphErrors* graphReq1oE       = new TGraphErrors(nP); 
  graphReq1oE->Sort();
  
  
  Int_t useEta[7][nEta]             = {{ 0}};
  Int_t useEtaE[7][nP]             = {{ 0}};
  Int_t usePhi[7][nPhi]             = {{ 0}};
  Int_t usePhiE[7][nP]             = {{ 0}};
  Int_t usePart[7]                  = { 0};
  Int_t usePartPhi[7]               = { 0};
  Int_t usePartEta[7]               = { 0};

  Int_t useParEnetPhi[7]            = { 0};
  Int_t useEnePhi[7][nP]         = {{ 0}};

  Int_t useParEnetEta[7]            = { 0};
  Int_t useEneEta[7][nP]         = {{ 0}};
          
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
  Int_t nActiveEtapartE[7]          = {0};
  Int_t nActivePhipart[7]           = {0};
  Int_t nActivePhipartE[7]          = {0};
  Int_t exampleBinsEtaE[4]          = {4, 6, 10, 15};
  Int_t exampleBinsPhiE[4]          = {4, 6, 10, 15};
  
  TFile* inputFile                      = new TFile(inputFileName.Data());
  TH3F* hist3DResoEta[7]                = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH2F* hist2DResoEta[7]                = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH2F* hist2DResoEtaVsE[7]             = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1F* hist1DResoEtabins3D[7][nP][nEta]= {{NULL}};
  TH1F* hist1DResoEtabins[7][nEta]      = {{NULL}};
  TH1F* hist1DResoEtaEbins[7][nP]       = {{NULL}};
  TH1F* histResoEtaVsEta[7]             = {NULL};
  TH1F* histMeanEtaVsEta[7]             = {NULL};
  TH1F* histSigmaEtaVsEta[7]            = {NULL};
  TH1F* histResoEtaVsE[7]               = {NULL};
  TH1F* histMeanEtaVsE[7]               = {NULL};
  TH1F* histSigmaEtaVsE[7]              = {NULL};
  TH1F* histResoEEtaVsEta[7][nP]        = {{NULL}};
  TH1F* histMeanEEtaVsEta[7][nP]        = {{NULL}};
  TH1F* histSigmaEEtaVsEta[7][nP]       = {{NULL}};
  TH1F* histResoEEtaVsE[7][nEta]        = {{NULL}};
  TH1F* histMeanEEtaVsE[7][nEta]        = {{NULL}};
  TH1F* histSigmaEEtaVsE[7][nEta]       = {{NULL}};

  TH3F* hist3DResoPhi[7]                  = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH2F* hist2DResoPhi[7]                  = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH2F* hist2DResoPhiVsE[7]               = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  TH1F* hist1DResoPhibins3D[7][nP][nPhi]  = {{NULL}};
  TH1F* hist1DResoPhibins[7][nPhi]        = {{NULL}};
  TH1F* hist1DResoPhiEbins[7][nP]         = {{NULL}};
  TH1F* histResoPhiVsPhi[7]               = {NULL};
  TH1F* histMeanPhiVsPhi[7]               = {NULL};
  TH1F* histResoPhiVsE[7]                 = {NULL};
  TH1F* histMeanPhiVsE[7]                 = {NULL};
  TH1F* histSigmaPhiVsE[7]                = {NULL};
  TH1F* histSigmaPhiVsPhi[7]              = {NULL};
  TH1F* histResoEPhiVsPhi[7][nP]          = {{NULL}};
  TH1F* histMeanEPhiVsPhi[7][nP]          = {{NULL}};
  TH1F* histSigmaEPhiVsPhi[7][nP]         = {{NULL}};
  TH1F* histResoEPhiVsE[7][nPhi]          = {{NULL}};
  TH1F* histMeanEPhiVsE[7][nPhi]          = {{NULL}};
  TH1F* histSigmaEPhiVsE[7][nPhi]         = {{NULL}};

  //========================================================================================================
  //================= Evaluate eta resolution ==============================================================
  //========================================================================================================  
  for (Int_t i = 0; i < nParticles; i++){
    TString tempName = Form("%s/h_CRH_EtaReso_%sE_Eta_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    cout << "reading: " << tempName << endl;
    hist3DResoEta[i]    = (TH3F*)inputFile->Get(tempName.Data());

    //*****************************************************************
    // create slices in E for eta res & evaluate
    //*****************************************************************
    if (hist3DResoEta[i]->GetEntries() > 10){
      hist3DResoEta[i]->GetYaxis()->SetRangeUser(minEtaMax,maxEtaMax);
      hist2DResoEtaVsE[i]    = (TH2F*)hist3DResoEta[i]->Project3D("zx");
      hist2DResoEtaVsE[i]->SetName(Form("projection_Eta_2D_%s_E", readNames[i].Data()));
      
      for (int ebin = 0; ebin < nP; ebin++){
        histResoEEtaVsEta[i][ebin]  = new TH1F(Form("h_%sreso_EEta_%i", readNames[i].Data(),ebin), "", nEta, partEtaCalo);
        histMeanEEtaVsEta[i][ebin]  = new TH1F(Form("h_%sreso_EEtaMean_%i", readNames[i].Data(),ebin), "", nEta, partEtaCalo);
        histSigmaEEtaVsEta[i][ebin] = new TH1F(Form("h_%sreso_EEtaSigma_%i", readNames[i].Data(),ebin), "", nEta, partEtaCalo);


        for (int etabin = minEtaBinFull[region]; etabin < maxEtaBinFull[region]+1; etabin++){
          if (!histResoEEtaVsE[i][etabin])histResoEEtaVsE[i][etabin]  = new TH1F(Form("h_%sreso_EEta_E_%i", readNames[i].Data(),etabin), "", nP, partP);
          if (!histMeanEEtaVsE[i][etabin])histMeanEEtaVsE[i][etabin]  = new TH1F(Form("h_%sreso_EEtaMean_E_%i", readNames[i].Data(),etabin), "", nP, partP);
          if (!histSigmaEEtaVsE[i][etabin])histSigmaEEtaVsE[i][etabin] = new TH1F(Form("h_%sreso_EEtaSigma_E_%i", readNames[i].Data(),etabin), "", nP, partP);
          
          hist1DResoEtabins3D[i][ebin][etabin]    = (TH1F*)hist3DResoEta[i]->ProjectionZ(Form("projection_Eta_3D_%s%d_%d", readNames[i].Data(), ebin, etabin), 
                                                                              hist3DResoEta[i]->GetXaxis()->FindBin(partP[ebin]), hist3DResoEta[i]->GetXaxis()->FindBin(partP[ebin+1]),
                                                                              hist3DResoEta[i]->GetYaxis()->FindBin(partEtaCalo[etabin]),hist3DResoEta[i]->GetYaxis()->FindBin(partEtaCalo[etabin+1]),"e");
          hist1DResoEtabins3D[i][ebin][etabin]->Rebin(rebinRes*4);
          if (hist1DResoEtabins3D[i][ebin][etabin]->GetEntries() > 10 ){
            hist1DResoEtabins3D[i][ebin][etabin]->Scale(1/hist1DResoEtabins3D[i][ebin][etabin]->GetMaximum());
            useParEnetEta[i]           = 1;
            useEneEta[i][ebin]          = 1;
            Double_t mean     = FindLargestBin1DHist(hist1DResoEtabins3D[i][ebin][etabin], -0.19, 0.19);
            Double_t meanErr  = hist1DResoEtabins3D[i][ebin][etabin]->GetMeanError();
            Double_t sigma    = hist1DResoEtabins3D[i][ebin][etabin]->GetRMS();
            Double_t sigmaErr = hist1DResoEtabins3D[i][ebin][etabin]->GetRMSError();

            Double_t reso     = sigma/mean;
            Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
            histResoEEtaVsEta[i][ebin]->SetBinContent(histResoEEtaVsEta[i][ebin]->FindBin(partEtaCalo[etabin]+0.01),  100*reso);
            histResoEEtaVsEta[i][ebin]->SetBinError(histResoEEtaVsEta[i][ebin]->FindBin(partEtaCalo[etabin]+0.01), resoErr);
            histMeanEEtaVsEta[i][ebin]->SetBinContent(histResoEEtaVsEta[i][ebin]->FindBin(partEtaCalo[etabin]+0.01), mean);
            histMeanEEtaVsEta[i][ebin]->SetBinError(histResoEEtaVsEta[i][ebin]->FindBin(partEtaCalo[etabin]+0.01), meanErr);
            histSigmaEEtaVsEta[i][ebin]->SetBinContent(histResoEEtaVsEta[i][ebin]->FindBin(partEtaCalo[etabin]+0.01), sigma);
            histSigmaEEtaVsEta[i][ebin]->SetBinError(histResoEEtaVsEta[i][ebin]->FindBin(partEtaCalo[etabin]+0.01), sigmaErr);  

            Int_t binE = histResoEEtaVsE[i][etabin]->FindBin(partP[ebin]+0.01);
            histResoEEtaVsE[i][etabin]->SetBinContent(binE,  100*reso);
            histResoEEtaVsE[i][etabin]->SetBinError(binE, resoErr);
            histMeanEEtaVsE[i][etabin]->SetBinContent(binE, mean);
            histMeanEEtaVsE[i][etabin]->SetBinError(binE, meanErr);
            histSigmaEEtaVsE[i][etabin]->SetBinContent(binE, sigma);
            histSigmaEEtaVsE[i][etabin]->SetBinError(binE, sigmaErr);  
          }      
        }
      }
    }
    //*****************************************************************
    // integrated eta res vs eta & E
    //*****************************************************************
    histResoEtaVsEta[i]       = new TH1F(Form("h_%sreso_Eta", readNames[i].Data()), "", nEta, partEtaCalo);
    histMeanEtaVsEta[i]       = new TH1F(Form("h_%sreso_EtaMean", readNames[i].Data()), "", nEta, partEtaCalo);
    histSigmaEtaVsEta[i]      = new TH1F(Form("h_%sreso_EtaSigma", readNames[i].Data()), "", nEta, partEtaCalo);
    histResoEtaVsE[i]         = new TH1F(Form("h_%sreso_Eta_E", readNames[i].Data()), "", nP, partP);
    histMeanEtaVsE[i]         = new TH1F(Form("h_%sreso_EtaMean_E", readNames[i].Data()), "", nP, partP);
    histSigmaEtaVsE[i]        = new TH1F(Form("h_%sreso_EtaSigma_E", readNames[i].Data()), "", nP, partP);

    tempName = Form("%s/h_CRH_EtaReso_%sEta_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    cout << "reading: " << tempName << endl;
    hist2DResoEta[i]    = (TH2F*)inputFile->Get(tempName.Data());

    if (hist2DResoEta[i]->GetEntries() > 0){
      for (int etabin = minEtaBinFull[region]; etabin < maxEtaBinFull[region]+1; etabin++){
        hist1DResoEtabins[i][etabin]    = (TH1F*)hist2DResoEta[i]->ProjectionY(Form("projection_dEta_%s%d", readNames[i].Data(), etabin),
                                                                              hist2DResoEta[i]->GetXaxis()->FindBin(partEtaCalo[etabin]), 
                                                                              hist2DResoEta[i]->GetXaxis()->FindBin(partEtaCalo[etabin+1]),"e");
        hist1DResoEtabins[i][etabin]->Rebin(rebinRes);
        if (hist1DResoEtabins[i][etabin]->GetEntries() > 10 ){
          hist1DResoEtabins[i][etabin]->Scale(1/hist1DResoEtabins[i][etabin]->GetMaximum());
          usePartEta[i]            = 1;
          useEta[i][etabin]         = 1;
          nActiveEtapart[i]++;
          
          Double_t mean     = FindLargestBin1DHist(hist1DResoEtabins[i][etabin], -0.19, 0.19);
          Double_t meanErr  = hist1DResoEtabins[i][etabin]->GetMeanError();
          Double_t sigma    = hist1DResoEtabins[i][etabin]->GetRMS();
          Double_t sigmaErr = hist1DResoEtabins[i][etabin]->GetRMSError();
          
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          histResoEtaVsEta[i]->SetBinContent(histResoEtaVsEta[i]->FindBin(partEtaCalo[etabin]+0.01),  100*reso);
          histResoEtaVsEta[i]->SetBinError(histResoEtaVsEta[i]->FindBin(partEtaCalo[etabin]+0.01), resoErr);
          histMeanEtaVsEta[i]->SetBinContent(histResoEtaVsEta[i]->FindBin(partEtaCalo[etabin]+0.01), mean);
          histMeanEtaVsEta[i]->SetBinError(histResoEtaVsEta[i]->FindBin(partEtaCalo[etabin]+0.01), meanErr);
          histSigmaEtaVsEta[i]->SetBinContent(histResoEtaVsEta[i]->FindBin(partEtaCalo[etabin]+0.01), sigma);
          histSigmaEtaVsEta[i]->SetBinError(histResoEtaVsEta[i]->FindBin(partEtaCalo[etabin]+0.01), sigmaErr);  
        }
      }
    }
    if (!hist2DResoEtaVsE[i]) continue;
    if (hist2DResoEtaVsE[i]->GetEntries() > 0){
      for (int ebin = 0; ebin < nP; ebin++){
        hist1DResoEtaEbins[i][ebin]    = (TH1F*)hist2DResoEtaVsE[i]->ProjectionY(Form("projection_dEta_E_%s%d", readNames[i].Data(), ebin),
                                                                              hist2DResoEtaVsE[i]->GetXaxis()->FindBin(partP[ebin]), 
                                                                              hist2DResoEtaVsE[i]->GetXaxis()->FindBin(partP[ebin+1]),"e");
        hist1DResoEtaEbins[i][ebin]->Rebin(rebinRes);
        if (hist1DResoEtaEbins[i][ebin]->GetEntries() > 10 ){
          hist1DResoEtaEbins[i][ebin]->Scale(1/hist1DResoEtaEbins[i][ebin]->GetMaximum());
          useEtaE[i][ebin]         = 1;
          nActiveEtapartE[i]++;
          
          Double_t mean     = FindLargestBin1DHist(hist1DResoEtaEbins[i][ebin], -0.19, 0.19);
          Double_t meanErr  = hist1DResoEtaEbins[i][ebin]->GetMeanError();
          Double_t sigma    = hist1DResoEtaEbins[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoEtaEbins[i][ebin]->GetRMSError();
          
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          histResoEtaVsE[i]->SetBinContent(histResoEtaVsE[i]->FindBin(partP[ebin]+0.01),  100*reso);
          histResoEtaVsE[i]->SetBinError(histResoEtaVsE[i]->FindBin(partP[ebin]+0.01), resoErr);
          histMeanEtaVsE[i]->SetBinContent(histResoEtaVsE[i]->FindBin(partP[ebin]+0.01), mean);
          histMeanEtaVsE[i]->SetBinError(histResoEtaVsE[i]->FindBin(partP[ebin]+0.01), meanErr);
          histSigmaEtaVsE[i]->SetBinContent(histResoEtaVsE[i]->FindBin(partP[ebin]+0.01), sigma);
          histSigmaEtaVsE[i]->SetBinError(histResoEtaVsE[i]->FindBin(partP[ebin]+0.01), sigmaErr);  
        }
      }
    }
  }
  //========================================================================================================
  //================= Evaluate phi resolution ==============================================================
  //========================================================================================================  
  for (Int_t i = 0; i < nParticles; i++){    
    TString tempName = Form("%s/h_CRH_PhiReso_%sE_Phi_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    cout << "reading: " << tempName << endl;
    hist3DResoPhi[i]    = (TH3F*)inputFile->Get(tempName.Data());

    //*****************************************************************
    // create slices in E for phi res & evaluate
    //*****************************************************************
    if (hist3DResoPhi[i]->GetEntries() > 10){
      hist3DResoPhi[i]->GetYaxis()->SetRangeUser(minPhiMax,maxPhiMax);
      hist2DResoPhiVsE[i]    = (TH2F*)hist3DResoPhi[i]->Project3D("zx");
      hist2DResoPhiVsE[i]->SetName(Form("projection_Phi_2D_%s_E", readNames[i].Data()));

      for (int ebin = 0; ebin < nP; ebin++){
        histResoEPhiVsPhi[i][ebin]  = new TH1F(Form("h_%sreso_EPhi_%i", readNames[i].Data(),ebin), "", nPhi, partPhi);
        histMeanEPhiVsPhi[i][ebin]  = new TH1F(Form("h_%sreso_EPhiMean_%i", readNames[i].Data(),ebin), "", nPhi, partPhi);
        histSigmaEPhiVsPhi[i][ebin] = new TH1F(Form("h_%sreso_EPhiSigma_%i", readNames[i].Data(),ebin), "", nPhi, partPhi);
        
        for (int phibin = 0; phibin < nPhi; phibin++){
          if (!histResoEPhiVsE[i][phibin])histResoEPhiVsE[i][phibin]  = new TH1F(Form("h_%sreso_EPhi_E_%i", readNames[i].Data(),phibin), "", nP, partP);
          if (!histMeanEPhiVsE[i][phibin])histMeanEPhiVsE[i][phibin]  = new TH1F(Form("h_%sreso_EPhiMean_E_%i", readNames[i].Data(),phibin), "", nP, partP);
          if (!histSigmaEPhiVsE[i][phibin])histSigmaEPhiVsE[i][phibin] = new TH1F(Form("h_%sreso_EPhiSigma_E_%i", readNames[i].Data(),phibin), "", nP, partP);

          hist1DResoPhibins3D[i][ebin][phibin]    = (TH1F*)hist3DResoPhi[i]->ProjectionZ(Form("projection_Phi_3D_%s%d_%d", readNames[i].Data(), ebin, phibin), 
                                                                              hist3DResoPhi[i]->GetXaxis()->FindBin(partP[ebin]), hist3DResoPhi[i]->GetXaxis()->FindBin(partP[ebin+1]),
                                                                              hist3DResoPhi[i]->GetYaxis()->FindBin(partPhi[phibin]),hist3DResoPhi[i]->GetYaxis()->FindBin(partPhi[phibin+1]),"e");
          hist1DResoPhibins3D[i][ebin][phibin]->Rebin(rebinRes*2);
          if (hist1DResoPhibins3D[i][ebin][phibin]->GetEntries() > 10 ){
            hist1DResoPhibins3D[i][ebin][phibin]->Scale(1/hist1DResoPhibins3D[i][ebin][phibin]->GetMaximum());
        
            useParEnetPhi[i]        = 1;
            useEnePhi[i][ebin]      = 1;
          
            Double_t mean     = FindLargestBin1DHist(hist1DResoPhibins3D[i][ebin][phibin], -0.8, 0.8);
            Double_t meanErr  = hist1DResoPhibins3D[i][ebin][phibin]->GetMeanError();
            Double_t sigma    = hist1DResoPhibins3D[i][ebin][phibin]->GetRMS();
            Double_t sigmaErr = hist1DResoPhibins3D[i][ebin][phibin]->GetRMSError();

            Double_t reso     = sigma/mean;
            Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
            histResoEPhiVsPhi[i][ebin]->SetBinContent(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[phibin]+0.001),  100*reso);
            histResoEPhiVsPhi[i][ebin]->SetBinError(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[phibin]+0.001), resoErr);
            histMeanEPhiVsPhi[i][ebin]->SetBinContent(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[phibin]+0.001), mean);
            histMeanEPhiVsPhi[i][ebin]->SetBinError(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[phibin]+0.001), meanErr);
            histSigmaEPhiVsPhi[i][ebin]->SetBinContent(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[phibin]+0.001), sigma);
            histSigmaEPhiVsPhi[i][ebin]->SetBinError(histResoEPhiVsPhi[i][ebin]->FindBin(partPhi[phibin]+0.001), sigmaErr);  
            
            Int_t binE = histResoEPhiVsE[i][phibin]->FindBin(partP[ebin]+0.01);
            histResoEPhiVsE[i][phibin]->SetBinContent(binE,  100*reso);
            histResoEPhiVsE[i][phibin]->SetBinError(binE, resoErr);
            histMeanEPhiVsE[i][phibin]->SetBinContent(binE, mean);
            histMeanEPhiVsE[i][phibin]->SetBinError(binE, meanErr);
            histSigmaEPhiVsE[i][phibin]->SetBinContent(binE, sigma);
            histSigmaEPhiVsE[i][phibin]->SetBinError(binE, sigmaErr);  
          }  
        }
      }
    }
    
    //*****************************************************************
    // integrated phi res vs phi & E
    //*****************************************************************
    histResoPhiVsPhi[i]       = new TH1F(Form("h_%sreso_Phi", readNames[i].Data()), "", nPhi, partPhi);
    histMeanPhiVsPhi[i]       = new TH1F(Form("h_%sreso_PhiMean", readNames[i].Data()), "", nPhi, partPhi);
    histSigmaPhiVsPhi[i]      = new TH1F(Form("h_%sreso_PhiSigma", readNames[i].Data()), "", nPhi, partPhi);
    histResoPhiVsE[i]         = new TH1F(Form("h_%sreso_Phi_E", readNames[i].Data()), "", nP, partP);
    histMeanPhiVsE[i]         = new TH1F(Form("h_%sreso_PhiMean_E", readNames[i].Data()), "", nP, partP);
    histSigmaPhiVsE[i]        = new TH1F(Form("h_%sreso_PhiSigma_E", readNames[i].Data()), "", nP, partP);

    TString tempNamePhi = Form("%s/h_CRH_PhiReso_%sPhi_%s_%s",caloNameRead.Data(), readNames[i].Data(), caloNameRead.Data(), clusterizerName.Data() );
    cout << "reading: " << tempNamePhi << endl;

    hist2DResoPhi[i]    = (TH2F*)inputFile->Get(tempNamePhi.Data());
    if (hist2DResoPhi[i]->GetEntries() > 0){
      for (int phibin = 0; phibin < nPhi; phibin++){
        hist1DResoPhibins[i][phibin]    = (TH1F*)hist2DResoPhi[i]->ProjectionY(Form("projection_dPhi_%s%d", readNames[i].Data(), phibin),  
                                                                               hist2DResoPhi[i]->GetXaxis()->FindBin(partPhi[phibin]), 
                                                                               hist2DResoPhi[i]->GetXaxis()->FindBin(partPhi[phibin+1]),"e");
        if (hist1DResoPhibins[i][phibin]->GetEntries() > 10 ){
          hist1DResoPhibins[i][phibin]->Scale(1/hist1DResoPhibins[i][phibin]->GetMaximum());
          usePartPhi[i]           = 1;
          usePhi[i][phibin]         = 1;
          nActivePhipart[i]++;
          Double_t mean     = FindLargestBin1DHist(hist1DResoPhibins[i][phibin], -0.8, 0.8);
          Double_t meanErr  = hist1DResoPhibins[i][phibin]->GetMeanError();
          Double_t sigma    = hist1DResoPhibins[i][phibin]->GetRMS();
          Double_t sigmaErr = hist1DResoPhibins[i][phibin]->GetRMSError();

          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          histResoPhiVsPhi[i]->SetBinContent(histResoPhiVsPhi[i]->FindBin(partPhi[phibin]+0.001),  100*reso);
          histResoPhiVsPhi[i]->SetBinError(histResoPhiVsPhi[i]->FindBin(partPhi[phibin]+0.001), resoErr);
          histMeanPhiVsPhi[i]->SetBinContent(histResoPhiVsPhi[i]->FindBin(partPhi[phibin]+0.001), mean);
          histMeanPhiVsPhi[i]->SetBinError(histResoPhiVsPhi[i]->FindBin(partPhi[phibin]+0.001), meanErr);
          histSigmaPhiVsPhi[i]->SetBinContent(histResoPhiVsPhi[i]->FindBin(partPhi[phibin]+0.001), sigma);
          histSigmaPhiVsPhi[i]->SetBinError(histResoPhiVsPhi[i]->FindBin(partPhi[phibin]+0.001), sigmaErr);
        }  
      }
    }
    if (!hist2DResoPhiVsE[i]) continue;
    if (hist2DResoPhiVsE[i]->GetEntries() > 0){
      for (int ebin = 0; ebin < nP; ebin++){
        hist1DResoPhiEbins[i][ebin]    = (TH1F*)hist2DResoPhiVsE[i]->ProjectionY(Form("projection_dPhi_E_%s%d", readNames[i].Data(), ebin),
                                                                              hist2DResoPhiVsE[i]->GetXaxis()->FindBin(partP[ebin]), 
                                                                              hist2DResoPhiVsE[i]->GetXaxis()->FindBin(partP[ebin+1]),"e");
        hist1DResoPhiEbins[i][ebin]->Rebin(rebinRes);
        if (hist1DResoPhiEbins[i][ebin]->GetEntries() > 10 ){
          hist1DResoPhiEbins[i][ebin]->Scale(1/hist1DResoPhiEbins[i][ebin]->GetMaximum());
          usePhiE[i][ebin]         = 1;
          nActivePhipartE[i]++;
          
          Double_t mean     = FindLargestBin1DHist(hist1DResoPhiEbins[i][ebin], -0.8, 0.8);
          Double_t meanErr  = hist1DResoPhiEbins[i][ebin]->GetMeanError();
          Double_t sigma    = hist1DResoPhiEbins[i][ebin]->GetRMS();
          Double_t sigmaErr = hist1DResoPhiEbins[i][ebin]->GetRMSError();
          
          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt(TMath::Power(sigmaErr/sigma,2)+TMath::Power(meanErr/mean,2));
          histResoPhiVsE[i]->SetBinContent(histResoPhiVsE[i]->FindBin(partP[ebin]+0.01),  100*reso);
          histResoPhiVsE[i]->SetBinError(histResoPhiVsE[i]->FindBin(partP[ebin]+0.01), resoErr);
          histMeanPhiVsE[i]->SetBinContent(histResoPhiVsE[i]->FindBin(partP[ebin]+0.01), mean);
          histMeanPhiVsE[i]->SetBinError(histResoPhiVsE[i]->FindBin(partP[ebin]+0.01), meanErr);
          histSigmaPhiVsE[i]->SetBinContent(histResoPhiVsE[i]->FindBin(partP[ebin]+0.01), sigma);
          histSigmaPhiVsE[i]->SetBinError(histResoPhiVsE[i]->FindBin(partP[ebin]+0.01), sigmaErr);  
        }
      }
    }
    
  }
  
  // ====================================================================================//
  // ======== Phi resolution drawing ====================================================//
  // ====================================================================================//
  for (Int_t i = 0; i < nParticles; i++){
    
    if (!usePartPhi[i]) continue;
    //************************************************************************************
    // Single reso bins phi
    //************************************************************************************
    TCanvas* cPNGPhi;
    split_canvas(cPNGPhi, Form("cPNGPhi%d",i), nActivePhipart[i]);
    Int_t padnum =0;
      
    for(Int_t phibin=0; phibin< nPhi; phibin++){
      if (!usePhi[i][phibin]) continue;
      cPNGPhi->cd(padnum+1);
      DrawVirtualPadSettings( cPNGPhi->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
    
      SetStyleHistoTH1ForGraphs(hist1DResoPhibins[i][phibin], "(#varphi^{rec} - #varphi^{MC})", "counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);

      hist1DResoPhibins[i][phibin]->GetXaxis()->SetRangeUser(-0.8,0.8);
      hist1DResoPhibins[i][phibin]->GetYaxis()->SetRangeUser(0.0,hist1DResoPhibins[i][phibin]->GetMaximum()*2);

      DrawGammaSetMarker(hist1DResoPhibins[i][phibin], 20, 1.5, kBlack, kBlack);      
      hist1DResoPhibins[i][phibin]->Draw("p,e");    
      
      if (padnum == 0){
        drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
       
      drawLatexAdd(Form("%0.2f < #varphi < %0.2f ",partPhi[phibin], partPhi[phibin+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
     
      Int_t bin = histMeanPhiVsPhi[i]->FindBin(partPhi[phibin]);
      cout << "phi: "<<  bin << "\t" << partPhi[phibin] << "\t" << histMeanPhiVsPhi[i]->GetBinContent(bin) << "\t" << histSigmaPhiVsPhi[i]->GetBinContent(bin) << endl;
      DrawGammaLines(histMeanPhiVsPhi[i]->GetBinContent(bin), 
                     histMeanPhiVsPhi[i]->GetBinContent(bin), 
                     0., 0.2*hist1DResoPhibins[i][phibin]->GetMaximum(), 1, kGray+2, 2);
      DrawGammaLines(histMeanPhiVsPhi[i]->GetBinContent(bin)-histSigmaPhiVsPhi[i]->GetBinContent(bin), 
                     histMeanPhiVsPhi[i]->GetBinContent(bin)-histSigmaPhiVsPhi[i]->GetBinContent(bin), 
                     0., 0.2*hist1DResoPhibins[i][phibin]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(histMeanPhiVsPhi[i]->GetBinContent(bin)+histSigmaPhiVsPhi[i]->GetBinContent(bin), 
                     histMeanPhiVsPhi[i]->GetBinContent(bin)+histSigmaPhiVsPhi[i]->GetBinContent(bin), 
                     0., 0.2*hist1DResoPhibins[i][phibin]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNGPhi->Print(Form("%s/%sphidist_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));
   
    //************************************************************************************
    // Single reso bins phi, differnent E bins on top
    //************************************************************************************
    padnum =0; 
    for(Int_t phibin=0; phibin< nPhi; phibin++){
      if (!usePhi[i][phibin]) continue;
      cPNGPhi->cd(padnum+1);
        for (Int_t ex = 0; ex < 4; ex++){
          DrawGammaSetMarker(hist1DResoPhibins3D[i][exampleBinsPhiE[ex]][phibin], markerStylePart[ex], 0.5*markerSizePart[ex], colorPart[ex], colorPart[ex]);    
          hist1DResoPhibins3D[i][exampleBinsPhiE[ex]][phibin]->Draw("PE,same");
        }
      padnum++;
    }
    cPNGPhi->Print(Form("%s/%sphidist_Edep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

    //************************************************************************************
    // Single reso bins phi in E bins for whole calo
    //************************************************************************************
    split_canvas(cPNGPhi, Form("cPNGPhi%d",i), nActivePhipartE[i]);
    padnum =0;
    for(Int_t ebin=0; ebin< nP; ebin++){
      if (!usePhiE[i][ebin]) continue;
      cPNGPhi->cd(padnum+1);
      DrawVirtualPadSettings( cPNGPhi->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
    
      SetStyleHistoTH1ForGraphs(hist1DResoPhiEbins[i][ebin], "(#varphi^{rec} - #varphi^{MC})", "counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
      hist1DResoPhiEbins[i][ebin]->GetXaxis()->SetRangeUser(-0.8,0.8);
      hist1DResoPhiEbins[i][ebin]->GetYaxis()->SetRangeUser(0.,hist1DResoPhiEbins[i][ebin]->GetMaximum()*2);        
      DrawGammaSetMarker(hist1DResoPhiEbins[i][ebin], 20, 1.5, kBlack, kBlack);      
      hist1DResoPhiEbins[i][ebin]->Draw("p,e");

      if (padnum == 0){
        drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.1f< #varphi< %1.1f",minPhiMax,maxPhiMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
      drawLatexAdd(Form("%0.1f < #it{E} (GeV) < %0.1f ",partP[ebin], partP[ebin+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
      Int_t bin = histMeanPhiVsE[i]->FindBin(partP[ebin]);
      cout << "E: "<<  bin << "\t" << partP[ebin] << "\t" << histMeanPhiVsE[i]->GetBinContent(bin) << "\t" << histSigmaPhiVsE[i]->GetBinContent(bin) << endl;
      DrawGammaLines(histMeanPhiVsE[i]->GetBinContent(bin), 
                    histMeanPhiVsE[i]->GetBinContent(bin), 
                    0., 0.2*hist1DResoPhiEbins[i][ebin]->GetMaximum(), 1, kGray+2, 2);
      DrawGammaLines(histMeanPhiVsE[i]->GetBinContent(bin)-histSigmaPhiVsE[i]->GetBinContent(bin), 
                    histMeanPhiVsE[i]->GetBinContent(bin)-histSigmaPhiVsE[i]->GetBinContent(bin), 
                    0., 0.2*hist1DResoPhiEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(histMeanPhiVsE[i]->GetBinContent(bin)+histSigmaPhiVsE[i]->GetBinContent(bin), 
                    histMeanPhiVsE[i]->GetBinContent(bin)+histSigmaPhiVsE[i]->GetBinContent(bin), 
                    0., 0.2*hist1DResoPhiEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNGPhi->Print(Form("%s/%sphidist_Ebins_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

    //************************************************************************************
    // Single reso bins phi, different E bin on top
    //************************************************************************************
    padnum= 0;
    for(Int_t ebin=0; ebin< nP; ebin++){
      if (!usePhiE[i][ebin]) continue;
      cPNGPhi->cd(padnum+1);
        for (int phibin = 0; phibin < nPhi; phibin++){
          if (! hist1DResoPhibins3D[i][ebin][phibin]) continue;
          cout << phibin << "\t" << hist1DResoPhibins3D[i][ebin][phibin]->GetMaximum() << "\t" << hist1DResoPhibins3D[i][ebin][phibin]->GetMean() << "\t" << hist1DResoPhibins3D[i][ebin][phibin]->GetRMS() << endl;;
          DrawGammaSetMarker(hist1DResoPhibins3D[i][ebin][phibin], markerStylePhi[phibin], 0.5*markerSizePhi[phibin], colorPhi[phibin], colorPhi[phibin]);    
          hist1DResoPhibins3D[i][ebin][phibin]->Draw("PE,same");
        }
      padnum++;
    }
    cPNGPhi->Print(Form("%s/%sphidist_Ebins_Phidep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));
   
    if (hist2DResoPhi[i]->GetEntries() > 0){

      //************************************************************************************
      // 2D plot res Phi vs Phi
      //************************************************************************************
      TCanvas* phidist = new TCanvas("phidist","",0,0,1100,1000);
      DrawGammaCanvasSettings( phidist, 0.095, 0.1, 0.01, 0.105);
        SetStyleHistoTH2ForGraphs(hist2DResoPhi[i], "#varphi^{MC}","(#varphi^{rec} - #varphi^{MC}) (rad)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
        hist2DResoPhi[i]->GetYaxis()->SetNdivisions(505,kTRUE);
        hist2DResoPhi[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DResoPhi[i]->Draw("colz");
      phidist->Print(Form("%s/%sphi_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

      //************************************************************************************
      // 2D plot res Phi vs E
      //************************************************************************************
      phidist->SetLogx();
        SetStyleHistoTH2ForGraphs(hist2DResoPhiVsE[i], "#it{E}^{MC} (GeV)","(#varphi^{rec} - #varphi^{MC}) (rad)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        hist2DResoPhiVsE[i]->GetYaxis()->SetNdivisions(505,kTRUE);
        hist2DResoPhiVsE[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DResoPhiVsE[i]->GetXaxis()->SetLabelOffset(-0.002);
        hist2DResoPhiVsE[i]->Draw("colz");
      phidist->Print(Form("%s/%sphi_E_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));


      //************************************************************************************
      // res Phi vs Phi mean integrated E
      //************************************************************************************
      TCanvas* phivar = new TCanvas("phivar","",0,0,1100,1000);
      DrawGammaCanvasSettings( phivar, 0.13, 0.01, 0.01, 0.105);
      phivar->cd();
        SetStyleHistoTH1ForGraphs(histMeanPhiVsPhi[i], "#varphi (rad)","#mu(#varphi^{rec} - #varphi^{MC}) (rad)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
        DrawGammaSetMarker(histMeanPhiVsPhi[i], 20, 1.5, kBlack, kBlack); 
        histMeanPhiVsPhi[i]->GetYaxis()->SetRangeUser(-0.14,0.14);
        histMeanPhiVsPhi[i]->Draw("PE");
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      phivar->SaveAs(Form("%s/%sphi-mean_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      //************************************************************************************
      // res Phi vs Phi mean different E
      //************************************************************************************      
      TLegend* legendEPhiMean  = GetAndSetLegend2(0.17, 0.95-(3*0.85*textSizeLabelsRel), 0.7, 0.95,0.85*textSizeLabelsPixel, 2, "", 43, 0.15);
      legendEPhiMean->AddEntry(histMeanPhiVsPhi[i],"#LT E #GT","p");
      histMeanPhiVsPhi[i]->GetYaxis()->SetRangeUser(-0.6,0.14);
      histMeanPhiVsPhi[i]->Draw("PE");
      drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

      for (Int_t ex = 0; ex < 4; ex++){
        DrawGammaSetMarker(histMeanEPhiVsPhi[i][exampleBinsPhiE[ex]], markerStylePart[ex], 0.5*markerSizePart[ex], colorPart[ex], colorPart[ex]);    
        histMeanEPhiVsPhi[i][exampleBinsPhiE[ex]]->Draw("PE,same");
        legendEPhiMean->AddEntry(histMeanEPhiVsPhi[i][exampleBinsPhiE[ex]], Form("%1.1f < #it{E} (GeV) < %1.1f", partP[exampleBinsPhiE[ex]], partP[exampleBinsPhiE[ex]+1]),"p");
      }
      legendEPhiMean->Draw();
      phivar->SaveAs(Form("%s/%sphi-mean_Edep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 
      
      //************************************************************************************
      // res Phi vs E mean different Phi
      //************************************************************************************
      TLegend* legendEPhiMeanE  = GetAndSetLegend2(0.17, 0.95-(nPhi/4*0.85*textSizeLabelsRel), 0.9, 0.95,0.65*textSizeLabelsPixel, 4, "", 43, 0.1);
      phivar->SetLogx();
      
        SetStyleHistoTH1ForGraphs(histMeanPhiVsE[i], "#it{E} (GeV)","#mu(#varphi^{rec} - #varphi^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.3, 510, 506);
        histMeanPhiVsE[i]->GetXaxis()->SetRangeUser(0.2,partP[nP]);
        histMeanPhiVsE[i]->GetYaxis()->SetRangeUser(-0.6,0.14);
        DrawGammaSetMarker(histMeanPhiVsE[i], 20, 1.5,  kBlack, kBlack);    
        histMeanPhiVsE[i]->Draw("PE");
        legendEPhiMeanE->AddEntry(histMeanPhiVsE[i], Form("%1.1f < #varphi < %1.1f", minPhiMax,maxPhiMax),"p");
        
        for (int phibin = 0; phibin < nPhi; phibin++){
          DrawGammaSetMarker(histMeanEPhiVsE[i][phibin], markerStylePhi[phibin], markerSizePhi[phibin], colorPhi[phibin], colorPhi[phibin]);    
          histMeanEPhiVsE[i][phibin]->Draw("PE,same");
          legendEPhiMeanE->AddEntry(histMeanEPhiVsE[i][phibin], Form("%1.1f < #varphi < %1.1f", partPhi[phibin],partPhi[phibin+1] ),"p");
        }
        legendEPhiMeanE->Draw();
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      phivar->SaveAs(Form("%s/%sphi-mean_E_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      
      //************************************************************************************
      // res Phi vs Phi sigma integrated E
      //************************************************************************************      
      phivar->cd();
      phivar->SetLogx(kFALSE);
        SetStyleHistoTH1ForGraphs(histSigmaPhiVsPhi[i], "#varphi (rad)","#sigma(#varphi^{rec} - #varphi^{MC}) (rad)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
        DrawGammaSetMarker(histSigmaPhiVsPhi[i], 20, 1.5, kBlack, kBlack);    
        histSigmaPhiVsPhi[i]->GetYaxis()->SetRangeUser(0,0.43);
        histSigmaPhiVsPhi[i]->Draw("PE");
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      phivar->SaveAs(Form("%s/%sphi-width_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      //************************************************************************************
      // res Phi vs Phi sigma different E
      //************************************************************************************            
      for (Int_t ex = 0; ex < 4; ex++){
        DrawGammaSetMarker(histSigmaEPhiVsPhi[i][exampleBinsPhiE[ex]], markerStylePart[ex], 0.5*markerSizePart[ex], colorPart[ex], colorPart[ex]);    
        histSigmaEPhiVsPhi[i][exampleBinsPhiE[ex]]->Draw("PE,same");
      }
      legendEPhiMean->Draw();
      phivar->SaveAs(Form("%s/%sphi-width_Edep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 
      
      //************************************************************************************
      // res Phi vs E sigma different Phi
      //************************************************************************************
      phivar->cd();
      phivar->SetLogx();
      
        SetStyleHistoTH1ForGraphs(histSigmaPhiVsE[i], "#it{E} (GeV)","#sigma(#varphi^{rec} - #varphi^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.3, 510, 506);
        histSigmaPhiVsE[i]->GetXaxis()->SetRangeUser(0.2,partP[nP]);
        histSigmaPhiVsE[i]->GetYaxis()->SetRangeUser(0,0.63);
        DrawGammaSetMarker(histSigmaPhiVsE[i], 20, 1.5,  kBlack, kBlack);    
        histSigmaPhiVsE[i]->Draw("PE");

        for (int phibin = 0; phibin < nPhi; phibin++){    
          DrawGammaSetMarker(histSigmaEPhiVsE[i][phibin], markerStylePhi[phibin], markerSizePhi[phibin], colorPhi[phibin], colorPhi[phibin]);    
          histSigmaEPhiVsE[i][phibin]->Draw("PE,same");
        }
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        legendEPhiMeanE->Draw();
      phivar->SaveAs(Form("%s/%sphi-width_E_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

    }

    // ====================================================================================//
    // ======== Eta resolution drawing ====================================================//
    // ====================================================================================//
    if (hist2DResoEta[i]->GetEntries() > 0){
      //************************************************************************************
      // Single reso bins eta
      //************************************************************************************
      if (!usePartEta[i]) continue;
      TCanvas* cPNG;
      split_canvas(cPNG, Form("cPNG%d",i), nActiveEtapart[i]);

      Int_t padnum =0;
      for(Int_t etabin=0; etabin< nEta; etabin++){

        if (!useEta[i][etabin]) continue;
        cPNG->cd(padnum+1);
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
      
        SetStyleHistoTH1ForGraphs(hist1DResoEtabins[i][etabin], "(#eta^{rec} - #eta^{MC})", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        hist1DResoEtabins[i][etabin]->GetXaxis()->SetRangeUser(-0.2,0.2);
        hist1DResoEtabins[i][etabin]->GetYaxis()->SetRangeUser(0,hist1DResoEtabins[i][etabin]->GetMaximum()*2);        
        DrawGammaSetMarker(hist1DResoEtabins[i][etabin], 20, 1.5, kBlack, kBlack);      
        hist1DResoEtabins[i][etabin]->Draw("p,e");

        if (padnum == 0){
          drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  //         drawLatexAdd(Form("%1.1f< #eta< %1.1f",minEtaMax,maxEtaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        drawLatexAdd(Form("%0.1f < #eta < %0.1f ",partEtaCalo[etabin], partEtaCalo[etabin+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      
        Int_t bin = histMeanEtaVsEta[i]->FindBin(partEtaCalo[etabin]);
        cout << "eta: "<<  bin << "\t" << partEtaCalo[etabin] << "\t" << histMeanEtaVsEta[i]->GetBinContent(bin) << "\t" << histSigmaEtaVsEta[i]->GetBinContent(bin) << endl;
        DrawGammaLines(histMeanEtaVsEta[i]->GetBinContent(bin), 
                      histMeanEtaVsEta[i]->GetBinContent(bin), 
                      0., 0.2*hist1DResoEtabins[i][etabin]->GetMaximum(), 1, kGray+2, 2);
        DrawGammaLines(histMeanEtaVsEta[i]->GetBinContent(bin)-histSigmaEtaVsEta[i]->GetBinContent(bin), 
                      histMeanEtaVsEta[i]->GetBinContent(bin)-histSigmaEtaVsEta[i]->GetBinContent(bin), 
                      0., 0.2*hist1DResoEtabins[i][etabin]->GetMaximum(), 1, kRed+2, 2);
        DrawGammaLines(histMeanEtaVsEta[i]->GetBinContent(bin)+histSigmaEtaVsEta[i]->GetBinContent(bin), 
                      histMeanEtaVsEta[i]->GetBinContent(bin)+histSigmaEtaVsEta[i]->GetBinContent(bin), 
                      0., 0.2*hist1DResoEtabins[i][etabin]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNG->Print(Form("%s/%setadist_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

      //************************************************************************************
      // Single reso bins eta, different E bin on top
      //************************************************************************************
      padnum= 0;
      for(Int_t etabin=0; etabin< nEta; etabin++){
        if (!useEta[i][etabin]) continue;
        cPNG->cd(padnum+1);
          for (Int_t ex = 0; ex < 4; ex++){
            DrawGammaSetMarker(hist1DResoEtabins3D[i][exampleBinsEtaE[ex]][etabin], markerStylePart[ex], 0.5*markerSizePart[ex], colorPart[ex], colorPart[ex]);    
            hist1DResoEtabins3D[i][exampleBinsEtaE[ex]][etabin]->Draw("PE,same");
          }
        padnum++;
      }
      cPNG->Print(Form("%s/%setadist_Edep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

      //************************************************************************************
      // Single reso bins eta in E bins for whole calo
      //************************************************************************************
      split_canvas(cPNG, Form("cPNG%d",i), nActiveEtapartE[i]);
      padnum =0;
      for(Int_t ebin=0; ebin< nP; ebin++){
        if (!useEtaE[i][ebin]) continue;
        cPNG->cd(padnum+1);
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
      
        SetStyleHistoTH1ForGraphs(hist1DResoEtaEbins[i][ebin], "(#eta^{rec} - #eta^{MC})", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        hist1DResoEtaEbins[i][ebin]->GetXaxis()->SetRangeUser(-0.2,0.2);
        hist1DResoEtaEbins[i][ebin]->GetYaxis()->SetRangeUser(0.,hist1DResoEtaEbins[i][ebin]->GetMaximum()*2);        
        DrawGammaSetMarker(hist1DResoEtaEbins[i][ebin], 20, 1.5, kBlack, kBlack);      
        hist1DResoEtaEbins[i][ebin]->Draw("p,e");

        if (padnum == 0){
          drawLatexAdd(Form("%s in %s", labelPart[i].Data() ,detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%1.1f< #eta< %1.1f",minEtaMax,maxEtaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        drawLatexAdd(Form("%0.1f < #it{E} (GeV) < %0.1f ",partP[ebin], partP[ebin+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      
        Int_t bin = histMeanEtaVsE[i]->FindBin(partP[ebin]);
        cout << "E: "<<  bin << "\t" << partP[ebin] << "\t" << histMeanEtaVsE[i]->GetBinContent(bin) << "\t" << histSigmaEtaVsE[i]->GetBinContent(bin) << endl;
        DrawGammaLines(histMeanEtaVsE[i]->GetBinContent(bin), 
                      histMeanEtaVsE[i]->GetBinContent(bin), 
                      0., 0.2*hist1DResoEtaEbins[i][ebin]->GetMaximum(), 1, kGray+2, 2);
        DrawGammaLines(histMeanEtaVsE[i]->GetBinContent(bin)-histSigmaEtaVsE[i]->GetBinContent(bin), 
                      histMeanEtaVsE[i]->GetBinContent(bin)-histSigmaEtaVsE[i]->GetBinContent(bin), 
                      0., 0.2*hist1DResoEtaEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);
        DrawGammaLines(histMeanEtaVsE[i]->GetBinContent(bin)+histSigmaEtaVsE[i]->GetBinContent(bin), 
                      histMeanEtaVsE[i]->GetBinContent(bin)+histSigmaEtaVsE[i]->GetBinContent(bin), 
                      0., 0.2*hist1DResoEtaEbins[i][ebin]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNG->Print(Form("%s/%setadist_Ebins_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

      //************************************************************************************
      // Single reso bins eta, different E bin on top
      //************************************************************************************
      padnum= 0;
      for(Int_t ebin=0; ebin< nP; ebin++){
        if (!useEtaE[i][ebin]) continue;
        cPNG->cd(padnum+1);
          for (int etabin = minEtaBinFull[region]; etabin < maxEtaBinFull[region]+1; etabin++){
            if (! hist1DResoEtabins3D[i][ebin][etabin]) continue;
            cout << etabin << "\t" << hist1DResoEtabins3D[i][ebin][etabin]->GetMaximum() << "\t" << hist1DResoEtabins3D[i][ebin][etabin]->GetMean() << "\t" << hist1DResoEtabins3D[i][ebin][etabin]->GetRMS() << endl;;
            DrawGammaSetMarker(hist1DResoEtabins3D[i][ebin][etabin], markerStyleEta[etabin], 0.5*markerSizeEta[etabin], colorEta[etabin], colorEta[etabin]);    
            hist1DResoEtabins3D[i][ebin][etabin]->Draw("PE,same");
          }
        padnum++;
      }
      cPNG->Print(Form("%s/%setadist_Ebins_Etadep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

      //************************************************************************************
      // 2D plot res Eta vs Eta
      //************************************************************************************
      TCanvas* etadist = new TCanvas("etadist","",0,0,1100,1000);
        DrawGammaCanvasSettings( etadist, 0.095, 0.1, 0.01, 0.105);
        SetStyleHistoTH2ForGraphs(hist2DResoEta[i], "#eta^{MC}","(#eta^{rec} - #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        hist2DResoEta[i]->GetYaxis()->SetNdivisions(505,kTRUE);
        hist2DResoEta[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DResoEta[i]->Draw("colz");
      etadist->Print(Form("%s/%seta_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

      //************************************************************************************
      // 2D plot res Eta vs E
      //************************************************************************************
      etadist->SetLogx();
        SetStyleHistoTH2ForGraphs(hist2DResoEtaVsE[i], "#it{E}^{MC} (GeV)","(#eta^{rec} - #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        hist2DResoEtaVsE[i]->GetYaxis()->SetNdivisions(505,kTRUE);
        hist2DResoEtaVsE[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DResoEtaVsE[i]->GetXaxis()->SetLabelOffset(-0.002);
        hist2DResoEtaVsE[i]->Draw("colz");
      etadist->Print(Form("%s/%seta_E_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data()));

      //************************************************************************************
      // res Eta vs Eta mean integrated E
      //************************************************************************************
      TCanvas* etavar = new TCanvas("etavar","",0,0,1100,1000);
      DrawGammaCanvasSettings( etavar, 0.13, 0.01, 0.01, 0.105);
      etavar->cd();
        SetStyleHistoTH1ForGraphs(histMeanEtaVsEta[i], "#eta","#mu(#eta^{rec} - #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3, 510, 508);
        histMeanEtaVsEta[i]->GetYaxis()->SetRangeUser(-0.1,0.1);
        DrawGammaSetMarker(histMeanEtaVsEta[i], 20, 1.5, kBlack, kBlack);    
        histMeanEtaVsEta[i]->Draw("PE");
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        
      etavar->SaveAs(Form("%s/%seta-mean_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      //************************************************************************************
      // res Eta vs Eta mean different E
      //************************************************************************************
      TLegend* legendEEtaMean  = GetAndSetLegend2(0.17, 0.95-(3*0.85*textSizeLabelsRel), 0.7, 0.95,0.85*textSizeLabelsPixel, 2, "", 43, 0.15);
      legendEEtaMean->AddEntry(histMeanEtaVsEta[i],"#LT E #GT","p");
      for (Int_t ex = 0; ex < 4; ex++){
        DrawGammaSetMarker(histMeanEEtaVsEta[i][exampleBinsEtaE[ex]], markerStylePart[ex], 0.5*markerSizePart[ex], colorPart[ex], colorPart[ex]);    
        histMeanEEtaVsEta[i][exampleBinsEtaE[ex]]->Draw("PE,same");
        legendEEtaMean->AddEntry(histMeanEEtaVsEta[i][exampleBinsEtaE[ex]], Form("%1.1f < #it{E} (GeV) < %1.1f", partP[exampleBinsEtaE[ex]], partP[exampleBinsEtaE[ex]+1]),"p");
      }
      legendEEtaMean->Draw();
      etavar->SaveAs(Form("%s/%seta-mean_Edep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      //************************************************************************************
      // res Eta vs E mean different Eta
      //************************************************************************************
      TLegend* legendEEtaMeanE  = GetAndSetLegend2(0.17, 0.95-(maxNEtaBinsFull[region]/2*0.85*textSizeLabelsRel), 0.7, 0.95,0.85*textSizeLabelsPixel, 2, "", 43, 0.15);
      etavar->SetLogx();
      
        SetStyleHistoTH1ForGraphs(histMeanEtaVsE[i], "#it{E} (GeV)","#mu(#eta^{rec} - #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.3, 510, 506);
        histMeanEtaVsE[i]->GetXaxis()->SetRangeUser(0.2,partP[nP]);
        histMeanEtaVsE[i]->GetYaxis()->SetRangeUser(-0.13,0.21);
        DrawGammaSetMarker(histMeanEtaVsE[i], 20, 1.5,  kBlack, kBlack);    
        histMeanEtaVsE[i]->Draw("PE");
        legendEEtaMeanE->AddEntry(histMeanEtaVsE[i], Form("%1.1f < #eta < %1.1f", minEtaMax,maxEtaMax),"p");
        
        for (int etabin = minEtaBinFull[region]; etabin < maxEtaBinFull[region]+1; etabin++){
          DrawGammaSetMarker(histMeanEEtaVsE[i][etabin], markerStyleEta[etabin], markerSizeEta[etabin], colorEta[etabin], colorEta[etabin]);    
          histMeanEEtaVsE[i][etabin]->Draw("PE,same");
          legendEEtaMeanE->AddEntry(histMeanEEtaVsE[i][etabin], Form("%1.1f < #eta < %1.1f", partEtaCalo[etabin],partEtaCalo[etabin+1] ),"p");
        }
        legendEEtaMeanE->Draw();
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      etavar->SaveAs(Form("%s/%seta-mean_E_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      //************************************************************************************
      // res Eta vs Eta sigma integrated E
      //************************************************************************************      
      etavar->cd();
      etavar->SetLogx(0);
        SetStyleHistoTH1ForGraphs(histSigmaEtaVsEta[i], "#eta","#sigma(#eta^{rec} - #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.3, 510, 506);
        histSigmaEtaVsEta[i]->GetYaxis()->SetRangeUser(0,0.18);
        DrawGammaSetMarker(histSigmaEtaVsEta[i], 20, 1.5, kBlack, kBlack);    
        histSigmaEtaVsEta[i]->Draw("PE");
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);        
      etavar->SaveAs(Form("%s/%seta-width_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      //************************************************************************************
      // res Eta vs Eta sigma different E
      //************************************************************************************      
      DrawGammaSetMarker(histSigmaEtaVsEta[i], 20, 1.5, kBlack, kBlack);    
      histSigmaEtaVsEta[i]->Draw("PE");
      drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);        
      for (Int_t ex = 0; ex < 4; ex++){
        DrawGammaSetMarker(histSigmaEEtaVsEta[i][exampleBinsEtaE[ex]], markerStylePart[ex], 0.5*markerSizePart[ex], colorPart[ex], colorPart[ex]);    
        histSigmaEEtaVsEta[i][exampleBinsEtaE[ex]]->Draw("PE,same");
      }
      legendEEtaMean->Draw();
      etavar->SaveAs(Form("%s/%seta-width_Edep_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 

      //************************************************************************************
      // res Eta vs E sigma different Eta
      //************************************************************************************
      etavar->cd();
      etavar->SetLogx();
      
        SetStyleHistoTH1ForGraphs(histSigmaEtaVsE[i], "#it{E} (GeV)","#sigma(#eta^{rec} - #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.3, 510, 506);
        histSigmaEtaVsE[i]->GetXaxis()->SetRangeUser(0.2,partP[nP]);
        histSigmaEtaVsE[i]->GetYaxis()->SetRangeUser(0,0.18);
        DrawGammaSetMarker(histSigmaEtaVsE[i], 20, 1.5,  kBlack, kBlack);    
        histSigmaEtaVsE[i]->Draw("PE");

        for (int etabin = minEtaBinFull[region]; etabin < maxEtaBinFull[region]+1; etabin++){    
          DrawGammaSetMarker(histSigmaEEtaVsE[i][etabin], markerStyleEta[etabin], markerSizeEta[etabin], colorEta[etabin], colorEta[etabin]);    
          histSigmaEEtaVsE[i][etabin]->Draw("PE,same");
        }
        drawLatexAdd(perfLabel,0.16,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%s in %s", labelPart[i].Data(),detLabel.Data()),0.16,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        legendEEtaMeanE->Draw();
      etavar->SaveAs(Form("%s/%seta-width_E_%s.%s", outputDir.Data(), readNames[i].Data(), caloNameRead.Data(), suffix.Data())); 
    }
  }  
}
