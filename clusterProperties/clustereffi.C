#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void NormalizeByBinWidth(TH1D * hist) {
//	hist->Scale(1,"width");
//	return;
  int nBins = hist->GetNbinsX();
  for (int i = 1; i <= nBins; i++) {
    double fBinWidth = hist->GetXaxis()->GetBinWidth(i);
    hist->SetBinContent(i,hist->GetBinContent(i)/fBinWidth);
    hist->SetBinError(i,hist->GetBinError(i)/fBinWidth);
  }
}

void clustereffi(
                            TString inputFileNameTracking   = "file.root",
                            TString inputFileNameCluster    = "file.root",
                            TString calo            = "FHCAL",
                            TString suffix          = "pdf",
                            TString addLabel        = "",
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
  if (trackCuts == 1){
    readTrackClass                  = "LI2";
    addLabel                        = addLabel+"-LI2";
    writeLabel                      = "LI2";
    labelPlotCuts                   = "#geq 2 tracker hits";
  } else if (trackCuts == 2){
    readTrackClass                  = "LI3";
    addLabel                        = addLabel+"-LI3";
    writeLabel                      = "LI3";
    labelPlotCuts                   = "#geq 3 tracker hits";
  }
    
  for (Int_t eR = 0; eR < 3; eR++) cout << "before: "<< labelEtaRange[eR].Data() << endl;
  
  TString collisionSystem = "Pythia 6, e+p, 10+250 GeV";
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addLabel.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }

  TString outputDir                 = Form("plots/%s/ClusterEffi%s",dateForOutput.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  for (Int_t iCl = 0; iCl < nClus; iCl++) gSystem->Exec("mkdir -p "+outputDir+"/"+calo+nameClus[iCl]);
  
  
  TString detLabel = GetCollisionEnergy(addLabel);
  Int_t nActiveEta            = 6;
  Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;

  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;
  Int_t nActiceCl     = 6;
  if (calo.CompareTo("FEMC") == 0){
    nActiceCl     = 8;
    enableClus[2] = 1;
    enableClus[4] = 1;
  }
  
  TH1D* h_spectra_MC_E[nPID][nEta+1]              = {{NULL}};
  TH1D* h_spectraTr_MC_E[nPID][nEta+1]            = {{NULL}};
  TH1D* h_spectraCl_rec_E[nPID][nEta+1][nClus]    = {{{NULL}}};
  TH1D* h_spectraCl_rec_MCE[nPID][nEta+1][nClus]  = {{{NULL}}};
  TH1D* h_spectraCl_recSE_E[nPID][nEta+1][nClus]    = {{{NULL}}};
  TH1D* h_spectraCl_recSE_MCE[nPID][nEta+1][nClus]  = {{{NULL}}};
  TH1D* h_spectraCl_matched_recSE_MCE[nPID][nEta+1][nClus]  = {{{NULL}}};
  
  TH1D* h_spectraReb_MC_E[nPID][nEta+1]           = {{NULL}};
  TH1D* h_spectraTrReb_MC_E[nPID][nEta+1]            = {{NULL}};
  TH1D* h_spectraClReb_rec_E[nPID][nEta+1][nClus]   = {{{NULL}}};
  TH1D* h_spectraClReb_rec_MCE[nPID][nEta+1][nClus] = {{{NULL}}};
  TH1D* h_spectraClReb_recSE_E[nPID][nEta+1][nClus]   = {{{NULL}}};
  TH1D* h_spectraClReb_recSE_MCE[nPID][nEta+1][nClus] = {{{NULL}}};
  TH1D* h_spectraClReb_matched_recSE_MCE[nPID][nEta+1][nClus] = {{{NULL}}};

  TH1D* h_effi_rec_E[nPID][nEta+1][nClus]       = {{{NULL}}};
  TH1D* h_effi_rec_MCE[nPID][nEta+1][nClus]     = {{{NULL}}};
  TH1D* h_effi_recSE_E[nPID][nEta+1][nClus]       = {{{NULL}}};
  TH1D* h_effi_recSE_MCE[nPID][nEta+1][nClus]     = {{{NULL}}};
  TH1D* h_TMeffi_recSE_MCE[nPID][nEta+1][nClus]     = {{{NULL}}};
  
  TH2F* h_trackMapMC_eta_E[nPID]                = {NULL};
  TH2F* h_trackMapTr_eta_E[nPID]                = {NULL};
  TH2F* h_clusterMapRec_eta_E[nPID][nClus]      = {{NULL}};
  TH2F* h_clusterMapMC_eta_E[nPID][nClus]       = {{NULL}};
  TH2F* h_clusterMapSERec_eta_E[nPID][nClus]    = {{NULL}};
  TH2F* h_clusterMapSEMC_eta_E[nPID][nClus]     = {{NULL}};
  TH2F* h_clusterMapSEMC_matched_eta_E[nPID][nClus]     = {{NULL}};
  
  TH1D* h_cluster_NTowerMean_E[nClus]           = {NULL};
  TH1D* h_cluster_NClMean_E[nClus]              = {NULL};
  TH1D* h_cluster_NClMean_MCE[nClus]              = {NULL};
  
  TFile* inputFileTR  = new TFile(inputFileNameTracking.Data());
  TFile* inputFileCL  = new TFile(inputFileNameCluster.Data());
  for (Int_t iCl = 0; iCl < nClus; iCl++){
    h_cluster_NTowerMean_E[iCl]             = (TH1D*)inputFileCL->Get(Form("h_CS_NTowerMean_%s_%s_E", calo.Data(), nameClus[iCl].Data()));
    h_cluster_NClMean_E[iCl]                = (TH1D*)inputFileCL->Get(Form("h_CS_NClMean_%s_%s_E", calo.Data(), nameClus[iCl].Data()));
    h_cluster_NClMean_MCE[iCl]              = (TH1D*)inputFileCL->Get(Form("h_CS_NClMean_%s_%s_MCE", calo.Data(), nameClus[iCl].Data()));
  }
  for (Int_t pid = 1; pid < nPID; pid++){
    h_trackMapMC_eta_E[pid]                  = (TH2F*)inputFileTR->Get(Form("h_%s_MC_E", partNameET[pid].Data()));
    h_trackMapMC_eta_E[pid]->Sumw2();
    h_trackMapTr_eta_E[pid]                  = (TH2F*)inputFileTR->Get(Form("h_%s_rec_trueE", partNameET[pid].Data()));
    h_trackMapTr_eta_E[pid]->Sumw2();
    
    for (Int_t iCl = 0; iCl < nClus; iCl++){
      cout << Form("h_clusterizer_clsspec_particle_E_eta_%s_%s_%s", calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()) << endl;
      h_clusterMapRec_eta_E[pid][iCl]         = (TH2F*)inputFileCL->Get(Form("h_clusterizer_clsspec_particle_E_eta_%s_%s_%s", calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
      h_clusterMapMC_eta_E[pid][iCl]          = (TH2F*)inputFileCL->Get(Form("h_clusterizer_clsspecMC_particle_E_eta_%s_%s_%s", calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
      h_clusterMapSERec_eta_E[pid][iCl]       = (TH2F*)inputFileCL->Get(Form("h_clusterizer_clsspecSE_particle_E_eta_%s_%s_%s", calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
      h_clusterMapSEMC_eta_E[pid][iCl]        = (TH2F*)inputFileCL->Get(Form("h_clusterizer_clsspecSEMC_particle_E_eta_%s_%s_%s", calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
      h_clusterMapSEMC_matched_eta_E[pid][iCl]= (TH2F*)inputFileCL->Get(Form("h_clusterizer_clsspecSEMC_matched_particle_E_eta_%s_%s_%s", calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
    
      h_clusterMapRec_eta_E[pid][iCl]->Sumw2();
      h_clusterMapMC_eta_E[pid][iCl]->Sumw2();
      h_clusterMapSERec_eta_E[pid][iCl]->Sumw2();
      h_clusterMapSEMC_eta_E[pid][iCl]->Sumw2();
      h_clusterMapSEMC_matched_eta_E[pid][iCl]->Sumw2();
    }
  }
  
  for (Int_t iEta=0; iEta<maxEtaBinFull[2]+1;iEta++){
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput)std::cout << Form("%1.1f<#eta<%1.1f",etaMin,etaMax)  << std::endl;
    
    for (Int_t pid = 1; pid < nPID; pid++){

      h_spectra_MC_E[pid][iEta]          = (TH1D*)h_trackMapMC_eta_E[pid]->ProjectionX(Form("spectraMC%s_MCE_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapMC_eta_E[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapMC_eta_E[pid]->GetYaxis()->FindBin(etaMax-0.001),"e");       
      h_spectraReb_MC_E[pid][iEta]       = (TH1D*)h_spectra_MC_E[pid][iEta]->Rebin(nP-2, Form("spectraRebMC%s_MCE_%d",partName[pid].Data(), iEta),
                                                                                          partP);
      NormalizeByBinWidth(h_spectra_MC_E[pid][iEta]);
      NormalizeByBinWidth(h_spectraReb_MC_E[pid][iEta]);
      DrawGammaSetMarker(h_spectra_MC_E[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraReb_MC_E[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      
      h_spectraTr_MC_E[pid][iEta]          = (TH1D*)h_trackMapTr_eta_E[pid]->ProjectionX(Form("spectraTRMC%s_MCE_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapTr_eta_E[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapTr_eta_E[pid]->GetYaxis()->FindBin(etaMax-0.001),"e");       
      h_spectraTrReb_MC_E[pid][iEta]       = (TH1D*)h_spectraTr_MC_E[pid][iEta]->Rebin(nP-2, Form("spectraTRRebMC%s_MCE_%d",partName[pid].Data(), iEta),
                                                                                          partP);
      NormalizeByBinWidth(h_spectraTr_MC_E[pid][iEta]);
      NormalizeByBinWidth(h_spectraTrReb_MC_E[pid][iEta]);
      DrawGammaSetMarker(h_spectraTr_MC_E[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraTrReb_MC_E[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      
      
      for (Int_t iCl = 0; iCl < nClus; iCl++){
        h_spectraCl_rec_E[pid][iEta][iCl]          = (TH1D*)h_clusterMapRec_eta_E[pid][iCl]->ProjectionX(Form("spectraClRec%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()), 
                                                                            h_clusterMapRec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapRec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
        h_spectraClReb_rec_E[pid][iEta][iCl]       = (TH1D*)h_spectraCl_rec_E[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecReb%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),
                                                                                            partP);
        NormalizeByBinWidth(h_spectraCl_rec_E[pid][iEta][iCl]);
        NormalizeByBinWidth(h_spectraClReb_rec_E[pid][iEta][iCl]);
        DrawGammaSetMarker(h_spectraCl_rec_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        DrawGammaSetMarker(h_spectraClReb_rec_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        
        h_spectraCl_rec_MCE[pid][iEta][iCl]          = (TH1D*)h_clusterMapMC_eta_E[pid][iCl]->ProjectionX(Form("spectraClRec%s_MCE_%d_%s",partName[pid].Data(), iEta,
                                                                                                                nameClus[iCl].Data()), 
                                                                            h_clusterMapRec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapRec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
        h_spectraClReb_rec_MCE[pid][iEta][iCl]       = (TH1D*)h_spectraCl_rec_MCE[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecReb%s_MCE_%d_%s",partName[pid].Data(), iEta, 
                                                                                                                    nameClus[iCl].Data()), partP);
        NormalizeByBinWidth(h_spectraCl_rec_E[pid][iEta][iCl]);
        NormalizeByBinWidth(h_spectraClReb_rec_MCE[pid][iEta][iCl]);
        DrawGammaSetMarker(h_spectraCl_rec_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        DrawGammaSetMarker(h_spectraCl_rec_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

        h_spectraCl_recSE_E[pid][iEta][iCl]          = (TH1D*)h_clusterMapSERec_eta_E[pid][iCl]->ProjectionX(Form("spectraClRecSE%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()), 
                                                                            h_clusterMapSERec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapSERec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
        h_spectraClReb_recSE_E[pid][iEta][iCl]       = (TH1D*)h_spectraCl_recSE_E[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecSEReb%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),
                                                                                            partP);
        NormalizeByBinWidth(h_spectraCl_recSE_E[pid][iEta][iCl]);
        NormalizeByBinWidth(h_spectraClReb_recSE_E[pid][iEta][iCl]);
        DrawGammaSetMarker(h_spectraCl_recSE_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        DrawGammaSetMarker(h_spectraClReb_recSE_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        
        h_spectraCl_recSE_MCE[pid][iEta][iCl]          = (TH1D*)h_clusterMapSEMC_eta_E[pid][iCl]->ProjectionX(Form("spectraClRecSE%s_MCE_%d_%s",partName[pid].Data(), iEta,
                                                                                                                nameClus[iCl].Data()), 
                                                                            h_clusterMapSEMC_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapSEMC_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
        h_spectraClReb_recSE_MCE[pid][iEta][iCl]       = (TH1D*)h_spectraCl_recSE_MCE[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecSEReb%s_MCE_%d_%s",partName[pid].Data(), iEta, 
                                                                                                                    nameClus[iCl].Data()), partP);
        NormalizeByBinWidth(h_spectraCl_recSE_E[pid][iEta][iCl]);
        NormalizeByBinWidth(h_spectraClReb_recSE_MCE[pid][iEta][iCl]);
        DrawGammaSetMarker(h_spectraCl_recSE_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        DrawGammaSetMarker(h_spectraCl_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

        h_spectraCl_matched_recSE_MCE[pid][iEta][iCl]          = (TH1D*)h_clusterMapSEMC_matched_eta_E[pid][iCl]->ProjectionX(Form("spectraClRecSE_matched%s_MCE_%d_%s",partName[pid].Data(), iEta,
                                                                                                                nameClus[iCl].Data()), 
                                                                            h_clusterMapSEMC_matched_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapSEMC_matched_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
        h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl]       = (TH1D*)h_spectraCl_matched_recSE_MCE[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecSEReb_matched%s_MCE_%d_%s",partName[pid].Data(), iEta, 
                                                                                                                    nameClus[iCl].Data()), partP);
        NormalizeByBinWidth(h_spectraCl_matched_recSE_MCE[pid][iEta][iCl]);
        NormalizeByBinWidth(h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl]);
        
        
        
        h_effi_rec_E[pid][iEta][iCl]             = (TH1D*)h_spectraClReb_rec_E[pid][iEta][iCl]->Clone(Form("effi%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
        h_effi_rec_E[pid][iEta][iCl]->Divide(h_spectraClReb_rec_E[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");
        h_effi_rec_MCE[pid][iEta][iCl]           = (TH1D*)h_spectraClReb_rec_MCE[pid][iEta][iCl]->Clone(Form("effi%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
        h_effi_rec_MCE[pid][iEta][iCl]->Divide(h_spectraClReb_rec_MCE[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");

        h_effi_recSE_E[pid][iEta][iCl]             = (TH1D*)h_spectraClReb_recSE_E[pid][iEta][iCl]->Clone(Form("effiSE%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
        h_effi_recSE_E[pid][iEta][iCl]->Divide(h_spectraClReb_recSE_E[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");
        h_effi_recSE_MCE[pid][iEta][iCl]           = (TH1D*)h_spectraClReb_recSE_MCE[pid][iEta][iCl]->Clone(Form("effiSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
        h_effi_recSE_MCE[pid][iEta][iCl]->Divide(h_spectraClReb_recSE_MCE[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");

        h_TMeffi_recSE_MCE[pid][iEta][iCl]           = (TH1D*)h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl]->Clone(Form("TMeffiSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
        h_TMeffi_recSE_MCE[pid][iEta][iCl]->Divide(h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl],h_spectraTrReb_MC_E[pid][iEta],1,1,"B");
      }
    }
  }
  
  // 1D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
  DrawGammaCanvasSettings( cReso, 0.085, 0.025, 0.02, 0.105);
  // cReso->SetLogz();
  cReso->SetLogx();
      
  TH2F* histoDummyEffiE   = new TH2F("histoDummyEffiE","histoDummyEffiE",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiE, "#it{E} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.87);
  histoDummyEffiE->GetXaxis()->SetNoExponent();
  histoDummyEffiE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiE->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.87);
  histoDummyEffiMCE->GetXaxis()->SetNoExponent();
  histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);

  TH2F* histoDummyEffiTMMCE   = new TH2F("histoDummyEffiTMMCE","histoDummyEffiTMMCE",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiTMMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.87);
  histoDummyEffiTMMCE->GetXaxis()->SetNoExponent();
  histoDummyEffiTMMCE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiTMMCE->GetXaxis()->SetMoreLogLabels(kTRUE);

  TLegend* legendEffiE      = GetAndSetLegend2(0.12, 0.94-(nActiveEta/2*0.75*textSizeLabelsRel), 0.4, 0.94,0.75*textSizeLabelsPixel, 2, "", 43, 0.2);
  TLegend* legendEffiPID   = GetAndSetLegend2(0.12,  0.94-(3*0.85*textSizeLabelsRel), 0.4, 0.94, 0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  TLegend* legendEffiCl   = GetAndSetLegend2(0.12,  0.94-(nActiceCl/2*0.85*textSizeLabelsRel), 0.4,  0.94, 0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  
  for (Int_t iCl = 0; iCl < nClus; iCl++){
    for (Int_t pid = 1; pid < nPID; pid++){
      histoDummyEffiE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      legendEffiE->Clear();
      for(Int_t iEta=minEtaBinFull[2]; iEta<maxEtaBinFull[2]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_rec_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_rec_E[pid][iEta][iCl]->Draw("same,p");
        if ( h_effi_rec_E[pid][iEta][iCl]->GetMaximum() > 0){
          if (iEta == nEta )
            legendEffiE->AddEntry(h_effi_rec_E[pid][iEta][iCl],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
          else 
            legendEffiE->AddEntry(h_effi_rec_E[pid][iEta][iCl],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
        }
      }
      legendEffiE->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s, %s clusterizer", partLabel[pid].Data(), calo.Data(), nameClus[iCl].Data() ),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s%s/Effi_E_%s.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(), suffix.Data()));

      histoDummyEffiE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      if (pid == 0 ) legendEffiE->Clear();
      for(Int_t iEta=minEtaBinFull[2]; iEta<maxEtaBinFull[2]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_recSE_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_recSE_E[pid][iEta][iCl]->Draw("same,p");
      }
      legendEffiE->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s, %s clusterizer, single entry", partLabel[pid].Data(), calo.Data(), nameClus[iCl].Data() ),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s%s/EffiSE_E_%s.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(),  partName[pid].Data(), suffix.Data()));

      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for(Int_t iEta=minEtaBinFull[2]; iEta<maxEtaBinFull[2]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_rec_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_rec_MCE[pid][iEta][iCl]->Draw("same,p");
      }
      legendEffiE->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s, %s clusterizer", partLabel[pid].Data(), calo.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      cReso->Print(Form("%s/%s%s/Effi_MCE_%s.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(),  suffix.Data()));

      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for(Int_t iEta=0; iEta<nEta+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");
      }
      legendEffiE->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s, %s clusterizer, single entry", partLabel[pid].Data(), calo.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      cReso->Print(Form("%s/%s%s/EffiSE_MCE_%s.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(),   suffix.Data()));

      histoDummyEffiTMMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for(Int_t iEta=0; iEta<nEta+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_TMeffi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");
      }
      legendEffiE->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s, %s clusterizer, single entry", partLabel[pid].Data(), calo.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      cReso->Print(Form("%s/%s%s/TMEffiSE_MCE_%s.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(),   suffix.Data()));

    }
    for(Int_t iEta=minEtaBinFull[2]; iEta<maxEtaBinFull[2]+1;iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }
      histoDummyEffiMCE->Draw();
      legendEffiPID->Clear();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t pid =1; pid < nPID; pid++){
        
        DrawGammaSetMarker(h_effi_rec_MCE[pid][iEta][iCl], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
        h_effi_rec_MCE[pid][iEta][iCl]->Draw("same,p");      
        legendEffiPID->AddEntry(h_effi_rec_MCE[pid][iEta][iCl],partLabel[pid].Data(),"p");
      }
      legendEffiPID->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s, %s clusterizer", etaMin, etaMax, calo.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s%s/EffiPID_MCE_%d_%d.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      
      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t pid =1; pid < nPID; pid++){
        
        DrawGammaSetMarker(h_effi_recSE_MCE[pid][iEta][iCl], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
        h_effi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
      }
      legendEffiPID->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s, %s clusterizer, single entry", etaMin, etaMax, calo.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s%s/EffiPIDSE_MCE_%d_%d.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

      histoDummyEffiTMMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t pid =1; pid < nPID; pid++){
        
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[pid][iEta][iCl], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
        h_TMeffi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
      }
      legendEffiPID->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s, %s clusterizer, single entry", etaMin, etaMax, calo.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s%s/TMEffiPIDSE_MCE_%d_%d.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    }
  }
  
  for (Int_t pid =1; pid < nPID; pid++){
    for (Int_t iEta=minEtaBinFull[2]; iEta<maxEtaBinFull[2]+1;iEta++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (iEta < nEta){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }
      histoDummyEffiMCE->Draw();
      legendEffiCl->Clear();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t iCl = 0; iCl < nClus; iCl++){
        if (!enableClus[iCl]) continue;
        DrawGammaSetMarker(h_effi_rec_MCE[pid][iEta][iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
        h_effi_rec_MCE[pid][iEta][iCl]->Draw("same,p");      
        legendEffiCl->AddEntry(h_effi_rec_MCE[pid][iEta][iCl],nameClus[iCl].Data(),"p");
      }
      legendEffiCl->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s in %s", etaMin, etaMax, partLabel[pid].Data(), calo.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/EffiClusterizer_MCE_%s_%s_%d_%d.%s", outputDir.Data(), calo.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      
      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t iCl = 0; iCl < nClus; iCl++){
        if (!enableClus[iCl]) continue;
        DrawGammaSetMarker(h_effi_recSE_MCE[pid][iEta][iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
        h_effi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
      }
      legendEffiCl->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s in %s, single entry", etaMin, etaMax, partLabel[pid].Data(), calo.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/EffiClusterizerSE_MCE_%s_%s_%d_%d.%s", outputDir.Data(), calo.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

      histoDummyEffiTMMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t iCl = 0; iCl < nClus; iCl++){
        if (!enableClus[iCl]) continue;
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[pid][iEta][iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
        h_TMeffi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
      }
      legendEffiCl->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s in %s, single entry", etaMin, etaMax, partLabel[pid].Data(), calo.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/TMEffiClusterizerSE_MCE_%s_%s_%d_%d.%s", outputDir.Data(), calo.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    }
  }
  
  TH2F* histoDummyNTowerMean   = new TH2F("histoDummyNTowerMean","histoDummyNTowerMean",1000,0, 200,1000,0.0, 40);
  SetStyleHistoTH2ForGraphs(histoDummyNTowerMean, "#it{E} (GeV)","#LT#it{N}_{tower}#GT", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNTowerMean->GetXaxis()->SetNoExponent();
  histoDummyNTowerMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNTowerMean->GetXaxis()->SetMoreLogLabels(kTRUE);

  TLegend* legendNTowerCl   = GetAndSetLegend2(0.12,  0.94-(nActiceCl/2*0.85*textSizeLabelsRel), 0.4,  0.94, 0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  if (calo.CompareTo("FHCAL") == 0) histoDummyNTowerMean->GetYaxis()->SetRangeUser(0, 30);
  histoDummyNTowerMean->Draw();
  for (Int_t iCl = 0; iCl < nClus; iCl++){
    if (!enableClus[iCl]) continue;
    DrawGammaSetMarker(h_cluster_NTowerMean_E[iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
    h_cluster_NTowerMean_E[iCl]->Draw("same,p");      
    legendNTowerCl->AddEntry(h_cluster_NTowerMean_E[iCl],nameClus[iCl].Data(),"p");
  }
  legendNTowerCl->Draw();
  drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("%s", calo.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/NTowerInCluster_%s_Mean_E.%s", outputDir.Data(), calo.Data(),  suffix.Data()));

  TH2F* histoDummyNClPerParMean   = new TH2F("histoDummyNClPerParMean","histoDummyNClPerParMean",1000,0, 200,1000,0.0, 15);
  SetStyleHistoTH2ForGraphs(histoDummyNClPerParMean, "#it{E} (GeV)","#LT#it{N}_{cl}/particle #GT", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNClPerParMean->GetXaxis()->SetNoExponent();
  histoDummyNClPerParMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNClPerParMean->GetXaxis()->SetMoreLogLabels(kTRUE);

  if (calo.CompareTo("FHCAL") == 0) histoDummyNClPerParMean->GetYaxis()->SetRangeUser(0, 10);
  histoDummyNClPerParMean->Draw();
  DrawGammaLines(0.2, 200, 1., 1., 2, kGray+2, 7);
  for (Int_t iCl = 0; iCl < nClus; iCl++){
    if (!enableClus[iCl]) continue;
    DrawGammaSetMarker(h_cluster_NClMean_E[iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
    h_cluster_NClMean_E[iCl]->Draw("same,p");      
  }
  legendNTowerCl->Draw();
  drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("%s", calo.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/NClusterPerParticle_%s_Mean_E.%s", outputDir.Data(), calo.Data(),  suffix.Data()));

  histoDummyNClPerParMean->SetXTitle("#it{E}_{MC} (GeV)");
  histoDummyNClPerParMean->Draw();
  DrawGammaLines(0.2, 200, 1., 1., 2, kGray+2, 7);
  for (Int_t iCl = 0; iCl < nClus; iCl++){
    if (!enableClus[iCl]) continue;
    DrawGammaSetMarker(h_cluster_NClMean_MCE[iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
    h_cluster_NClMean_MCE[iCl]->Draw("same,p");      
  }
  legendNTowerCl->Draw();
  drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("%s", calo.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/NClusterPerParticle_%s_Mean_MCE.%s", outputDir.Data(), calo.Data(),  suffix.Data()));
  
//    
  TFile* outputFile  = new TFile(inputFileNameCluster.Data(),"UPDATE");
  TDirectoryFile* directoryClEffi = (TDirectoryFile*)outputFile->Get(Form("ClusterEfficiencyProjections_%s", calo.Data()));
  if (!directoryClEffi){
    outputFile->mkdir(Form("ClusterEfficiencyProjections_%s", calo.Data()));
    directoryClEffi = (TDirectoryFile*)outputFile->Get(Form("ClusterEfficiencyProjections_%s", calo.Data()));
  }
  directoryClEffi->cd();
  for (Int_t iCl = 0; iCl < nClus; iCl++){
    for (Int_t iEta=minEtaBinFull[2]; iEta<maxEtaBinFull[2]+1;iEta++){
      for (Int_t pid = 1; pid < 6; pid++){
        h_effi_rec_E[pid][iEta][iCl]->Write(Form("effi%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        h_effi_rec_MCE[pid][iEta][iCl]->Write(Form("effi%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        h_effi_recSE_E[pid][iEta][iCl]->Write(Form("effiSE%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        h_effi_recSE_MCE[pid][iEta][iCl]->Write(Form("effiSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
      }
    }
    h_cluster_NTowerMean_E[iCl]->Write(Form("h_CS_NTowerMean_%s_E", nameClus[iCl].Data()),TObject::kOverwrite);
    h_cluster_NClMean_E[iCl]->Write(Form("h_CS_NClMean_%s_E", nameClus[iCl].Data()),TObject::kOverwrite);
    h_cluster_NClMean_MCE[iCl]->Write(Form("h_CS_NClMean_%s_MCE", nameClus[iCl].Data()),TObject::kOverwrite);
  }
  outputFile->Write();
  outputFile->Close();
}
