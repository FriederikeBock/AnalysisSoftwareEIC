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
                            TString collisionsSys   = "PythiaMB",
                            TString addLabel        = "",
                            Int_t trackCuts         = 0,
                            unsigned short primaryTrackSource = 0
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
  
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addLabel.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }
  
  TString labelEnergy = "";
  if (collisionsSys.CompareTo("PythiaMB") == 0){
    labelEnergy   = "e-p: 18#times 275 GeV";
  } else if (collisionsSys.CompareTo("SingleElectron") == 0){
    labelEnergy   = "e^{-}";
  } else if (collisionsSys.CompareTo("SingleKaon") == 0){
    labelEnergy   = "K^{-}";
  } else if (collisionsSys.CompareTo("SinglePion") == 0){
    labelEnergy   = "#pi^{-}";
  } else if (collisionsSys.CompareTo("SingleProton") == 0){
    labelEnergy   = "p";
  } else if (collisionsSys.CompareTo("SingleNeutron") == 0){
    labelEnergy   = "n";
  } else if (collisionsSys.CompareTo("SinglePart") == 0){
    // labelEnergy   = "ECCE, single particle simulation";
    labelEnergy   = "#it{#bf{ECCE}} simulation";
  }
  
  Bool_t enableRecE = 0;
  if (!collisionsSys.Contains("Single")){
    enableRecE = 1;
  }
  Int_t nClusProcess = 3;
  
  TString outputDir                 = Form("plots/%s/ClusterEffi%s",dateForOutput.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  
  
  TString detLabel = GetCollisionEnergy(addLabel);
  Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;

  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;
  
  Int_t nActiceCl     = 6;
  Int_t nMaxClPart    = 10;
  Int_t nMaxTowCl     = 30;
  Int_t region        = 2;
  TString caloPlot = calo.Data();
  if (calo.CompareTo("FEMC") == 0){
    nActiceCl     = 2;
    enableClus[2] = 1;
    enableClus[4] = 1;
    nMaxClPart    = 5;
    nMaxTowCl     = 80;
    region        = 2;
  } else if (calo.CompareTo("LFHCAL") == 0){
    nActiceCl     = 2;
    enableClus[2] = 1;
    enableClus[4] = 1;
    nMaxClPart    = 5;
    nMaxTowCl     = 250;
    region        = 2;
  } else if (calo.CompareTo("BECAL") == 0){
    nActiceCl     = 2;
    enableClus[2] = 1;
    enableClus[4] = 1;
    nMaxClPart    = 5;
    nMaxTowCl     = 70;
    region        = 1;
    caloPlot      = "BEMC";
  } else if (calo.CompareTo("HCALIN") == 0){
    nActiceCl     = 2;
    enableClus[2] = 1;
    enableClus[4] = 1;
    nMaxClPart    = 5;
    nMaxTowCl     = 8;
    region        = 1;
    caloPlot      = "IHCAL";
  } else if (calo.CompareTo("HCALOUT") == 0){
    nActiceCl     = 2;
    enableClus[2] = 1;
    enableClus[4] = 1;
    nMaxClPart    = 5;
    nMaxTowCl     = 20;
    region        = 1;
    caloPlot      = "OHCAL";
  } else if (calo.CompareTo("EEMC") == 0){
    nActiceCl     = 2;
    enableClus[2] = 1;
    enableClus[4] = 1;
    nMaxClPart    = 5;
    nMaxTowCl     = 70;
    region        = 0;
  } else if (calo.CompareTo("EHCAL") == 0){
    nActiceCl     = 2;
    enableClus[2] = 1;
    enableClus[4] = 1;
    nMaxClPart    = 5;
    nMaxTowCl     = 70;
    region        = 0;
  }
  for (Int_t iCl = 0; iCl < nClusProcess; iCl++) gSystem->Exec("mkdir -p "+outputDir+"/"+caloPlot+nameClus[iCl]);

  Int_t nActiveEta            = maxNEtaBinsFull[region]+1;

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
  TH1D* h_TMeffiCls_recSE_MCE[nPID][nEta+1][nClus]     = {{{NULL}}};
  
  TH2F* h_trackMapMC_eta_E[nPID]                = {NULL};
  TH2F* h_trackMapTr_eta_E[nPID]                = {NULL};
  TH2F* h_clusterMapRec_eta_E[nPID][nClus]      = {{NULL}};
  TH2F* h_clusterMapMC_eta_E[nPID][nClus]       = {{NULL}};
  TH2F* h_clusterMapSERec_eta_E[nPID][nClus]    = {{NULL}};
  TH2F* h_clusterMapSEMC_eta_E[nPID][nClus]     = {{NULL}};
  TH2F* h_clusterMapSEMC_matched_eta_E[nPID][nClus]     = {{NULL}};
  
  TH1D* h_cluster_NTowerMean_E[nClus]           = {NULL};
  TH1D* h_cluster_NClMean_E[nClus]              = {NULL};
  TH1D* h_cluster_NClMean_MCE[nClus]            = {NULL};
  
  Bool_t enableParticle[nPID]                   = {0};
  Bool_t enableTM                               = 0;
  
  TFile* inputFileTR  = new TFile(inputFileNameTracking.Data());
  TFile* inputFileCL  = new TFile(inputFileNameCluster.Data());
  for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
    h_cluster_NTowerMean_E[iCl]             = (TH1D*)inputFileCL->Get(Form("%s/h_CS_NTowerMean_%s_%s_E", calo.Data(), calo.Data(), nameClus[iCl].Data()));
    h_cluster_NClMean_E[iCl]                = (TH1D*)inputFileCL->Get(Form("%s/h_CS_NClMean_%s_%s_E", calo.Data(), calo.Data(), nameClus[iCl].Data()));
    h_cluster_NClMean_MCE[iCl]              = (TH1D*)inputFileCL->Get(Form("%s/h_CS_NClMean_%s_%s_MCE", calo.Data(), calo.Data(), nameClus[iCl].Data()));
  }
  for (Int_t pid = 1; pid < nPID; pid++){
    h_trackMapMC_eta_E[pid]                   = (TH2F*)inputFileTR->Get(Form("h_%s_MC_E", partNameET[pid].Data()));
    h_trackMapMC_eta_E[pid]->Sumw2();
    if (h_trackMapMC_eta_E[pid]->GetEntries() > 0)
      enableParticle[pid]                    = 1;
    
    if (enableParticle[pid]){
      h_trackMapTr_eta_E[pid]                   = (TH2F*)inputFileTR->Get(Form("h_%s_rec_trueE_%u", partNameET[pid].Data(), primaryTrackSource));
      h_trackMapTr_eta_E[pid]->Sumw2();
      
      for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
        cout << Form("%s/h_clusterizer_clsspec_particle_E_eta_%s_%s_%s", calo.Data(), calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()) << endl;
        if (enableRecE){
          h_clusterMapRec_eta_E[pid][iCl]         = (TH2F*)inputFileCL->Get(Form("%s/h_clusterizer_clsspec_particle_E_eta_%s_%s_%s", calo.Data(), calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
          h_clusterMapSERec_eta_E[pid][iCl]       = (TH2F*)inputFileCL->Get(Form("%s/h_clusterizer_clsspecSE_particle_E_eta_%s_%s_%s", calo.Data(), calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
          h_clusterMapRec_eta_E[pid][iCl]->Sumw2();
          h_clusterMapSERec_eta_E[pid][iCl]->Sumw2();
        }
        h_clusterMapMC_eta_E[pid][iCl]          = (TH2F*)inputFileCL->Get(Form("%s/h_clusterizer_clsspecMC_particle_E_eta_%s_%s_%s", calo.Data(), calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
        if (h_trackMapTr_eta_E[pid]->GetEntries() == 0 && h_clusterMapMC_eta_E[pid][iCl]->GetEntries() == 0 ) // no rec particles -> probably only secondaries generated
          enableParticle[pid]                   = 0;
        h_clusterMapSEMC_eta_E[pid][iCl]        = (TH2F*)inputFileCL->Get(Form("%s/h_clusterizer_clsspecSEMC_particle_E_eta_%s_%s_%s", calo.Data(), calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
        h_clusterMapSEMC_matched_eta_E[pid][iCl]= (TH2F*)inputFileCL->Get(Form("%s/h_clusterizer_clsspecSEMC_matched_particle_E_eta_%s_%s_%s", calo.Data(), calo.Data(), nameClus[iCl].Data(), partNameET[pid].Data()));
        h_clusterMapMC_eta_E[pid][iCl]->Sumw2();
        h_clusterMapSEMC_eta_E[pid][iCl]->Sumw2();        
        h_clusterMapSEMC_matched_eta_E[pid][iCl]->Sumw2();
        
        if (h_clusterMapSEMC_matched_eta_E[pid][iCl]->GetEntries() > 0)
          enableTM  = 1;
        
      }
    }
  }
  
  for (Int_t iEta=0; iEta<maxEtaBinCaloDis[region]+1;iEta++){
    Double_t etaMin = partEta[minEtaBinCaloDis[region]];
    Double_t etaMax = partEta[maxEtaBinCaloDis[region]];
    if (iEta < maxEtaBinCaloDis[region]){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput)std::cout << Form("%1.1f < #it{#eta} < %1.1f",etaMin,etaMax)  << std::endl;
    
    for (Int_t pid = 1; pid < nPID; pid++){
      if (!enableParticle[pid]) continue;
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
      
      for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
        if (enableRecE){
          h_spectraCl_rec_E[pid][iEta][iCl]          = (TH1D*)h_clusterMapRec_eta_E[pid][iCl]->ProjectionX(Form("spectraClRec%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()), 
                                                                              h_clusterMapRec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapRec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
          h_spectraClReb_rec_E[pid][iEta][iCl]       = (TH1D*)h_spectraCl_rec_E[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecReb%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),
                                                                                              partP);
          NormalizeByBinWidth(h_spectraCl_rec_E[pid][iEta][iCl]);
          NormalizeByBinWidth(h_spectraClReb_rec_E[pid][iEta][iCl]);
          DrawGammaSetMarker(h_spectraCl_rec_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
          DrawGammaSetMarker(h_spectraClReb_rec_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        
          h_spectraCl_recSE_E[pid][iEta][iCl]          = (TH1D*)h_clusterMapSERec_eta_E[pid][iCl]->ProjectionX(Form("spectraClRecSE%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()), 
                                                                              h_clusterMapSERec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapSERec_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
          h_spectraClReb_recSE_E[pid][iEta][iCl]       = (TH1D*)h_spectraCl_recSE_E[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecSEReb%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),
                                                                                              partP);
          NormalizeByBinWidth(h_spectraCl_recSE_E[pid][iEta][iCl]);
          NormalizeByBinWidth(h_spectraClReb_recSE_E[pid][iEta][iCl]);
          DrawGammaSetMarker(h_spectraCl_recSE_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
          DrawGammaSetMarker(h_spectraClReb_recSE_E[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

          h_effi_rec_E[pid][iEta][iCl]             = (TH1D*)h_spectraClReb_rec_E[pid][iEta][iCl]->Clone(Form("effi%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
          h_effi_rec_E[pid][iEta][iCl]->Divide(h_spectraClReb_rec_E[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");
          h_effi_recSE_E[pid][iEta][iCl]             = (TH1D*)h_spectraClReb_recSE_E[pid][iEta][iCl]->Clone(Form("effiSE%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
          h_effi_recSE_E[pid][iEta][iCl]->Divide(h_spectraClReb_recSE_E[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");
        }
        h_spectraCl_rec_MCE[pid][iEta][iCl]          = (TH1D*)h_clusterMapMC_eta_E[pid][iCl]->ProjectionX(Form("spectraClRec%s_MCE_%d_%s",partName[pid].Data(), iEta,
                                                                                                                nameClus[iCl].Data()), 
                                                                            h_clusterMapMC_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapMC_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
        h_spectraClReb_rec_MCE[pid][iEta][iCl]       = (TH1D*)h_spectraCl_rec_MCE[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecReb%s_MCE_%d_%s",partName[pid].Data(), iEta, 
                                                                                                                    nameClus[iCl].Data()), partP);
        NormalizeByBinWidth(h_spectraCl_rec_MCE[pid][iEta][iCl]);
        DrawGammaSetMarker(h_spectraCl_rec_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        NormalizeByBinWidth(h_spectraClReb_rec_MCE[pid][iEta][iCl]);
        DrawGammaSetMarker(h_spectraClReb_rec_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

        h_spectraCl_recSE_MCE[pid][iEta][iCl]          = (TH1D*)h_clusterMapSEMC_eta_E[pid][iCl]->ProjectionX(Form("spectraClRecSE%s_MCE_%d_%s",partName[pid].Data(), iEta,
                                                                                                                nameClus[iCl].Data()), 
                                                                            h_clusterMapSEMC_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapSEMC_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
        h_spectraClReb_recSE_MCE[pid][iEta][iCl]       = (TH1D*)h_spectraCl_recSE_MCE[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecSEReb%s_MCE_%d_%s",partName[pid].Data(), iEta, 
                                                                                                                    nameClus[iCl].Data()), partP);
        NormalizeByBinWidth(h_spectraCl_recSE_MCE[pid][iEta][iCl]);
        NormalizeByBinWidth(h_spectraClReb_recSE_MCE[pid][iEta][iCl]);
        DrawGammaSetMarker(h_spectraCl_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        DrawGammaSetMarker(h_spectraClReb_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
       
        h_effi_rec_MCE[pid][iEta][iCl]           = (TH1D*)h_spectraClReb_rec_MCE[pid][iEta][iCl]->Clone(Form("effi%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
        h_effi_rec_MCE[pid][iEta][iCl]->Divide(h_spectraClReb_rec_MCE[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");

        h_effi_recSE_MCE[pid][iEta][iCl]           = (TH1D*)h_spectraClReb_recSE_MCE[pid][iEta][iCl]->Clone(Form("effiSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
        h_effi_recSE_MCE[pid][iEta][iCl]->Divide(h_spectraClReb_recSE_MCE[pid][iEta][iCl],h_spectraReb_MC_E[pid][iEta],1,1,"B");

        if (enableTM){
          h_spectraCl_matched_recSE_MCE[pid][iEta][iCl]          = (TH1D*)h_clusterMapSEMC_matched_eta_E[pid][iCl]->ProjectionX(Form("spectraClRecSE_matched%s_MCE_%d_%s",partName[pid].Data(), iEta,
                                                                                                                  nameClus[iCl].Data()), 
                                                                              h_clusterMapSEMC_matched_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMin+0.001), h_clusterMapSEMC_matched_eta_E[pid][iCl]->GetYaxis()->FindBin(etaMax-0.001),"e");       
          h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl]       = (TH1D*)h_spectraCl_matched_recSE_MCE[pid][iEta][iCl]->Rebin(nP-2, Form("spectraClRecSEReb_matched%s_MCE_%d_%s",partName[pid].Data(), iEta, 
                                                                                                                      nameClus[iCl].Data()), partP);
          NormalizeByBinWidth(h_spectraCl_matched_recSE_MCE[pid][iEta][iCl]);
          NormalizeByBinWidth(h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl]);

          h_TMeffi_recSE_MCE[pid][iEta][iCl]           = (TH1D*)h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl]->Clone(Form("TMeffiSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
          h_TMeffi_recSE_MCE[pid][iEta][iCl]->Divide(h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl],h_spectraTrReb_MC_E[pid][iEta],1,1,"B");
          h_TMeffiCls_recSE_MCE[pid][iEta][iCl]           = (TH1D*)h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl]->Clone(Form("TMeffiClsSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()));
          h_TMeffiCls_recSE_MCE[pid][iEta][iCl]->Divide(h_spectraClReb_matched_recSE_MCE[pid][iEta][iCl],h_spectraClReb_recSE_MCE[pid][iEta][iCl],1,1,"B");
        }
      }
    }
  }
  
  // 1D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
  DrawGammaCanvasSettings( cReso, 0.085, 0.025, 0.02, 0.105);
  // cReso->SetLogz();
  cReso->SetLogx();
      
  TH2F* histoDummyEffiE   = new TH2F("histoDummyEffiE","histoDummyEffiE",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiE, "#it{E} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.18*textSizeSinglePad, 0.9,0.7);
  histoDummyEffiE->GetXaxis()->SetNoExponent();
  histoDummyEffiE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiE->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.18*textSizeSinglePad, 0.9,0.7);
  histoDummyEffiMCE->GetXaxis()->SetNoExponent();
  histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);

  TH2F* histoDummyEffiTMMCE   = new TH2F("histoDummyEffiTMMCE","histoDummyEffiTMMCE",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiTMMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.18*textSizeSinglePad, 0.9,0.7);
  histoDummyEffiTMMCE->GetXaxis()->SetNoExponent();
  histoDummyEffiTMMCE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiTMMCE->GetXaxis()->SetMoreLogLabels(kTRUE);

  TLegend* legendEffiE      = GetAndSetLegend2(0.12, 0.94-(nActiveEta/2*0.85*textSizeLabelsRel), 0.5, 0.94,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
  TLegend* legendEffiPID    = GetAndSetLegend2(0.12,  0.94-(3*0.85*textSizeLabelsRel), 0.5, 0.94, 0.85*textSizeLabelsRel, 2, "", 42, 0.25);
  TLegend* legendEffiCl     = GetAndSetLegend2(0.12,  0.94-(nActiceCl/2*0.85*textSizeLabelsRel), 0.5,  0.94, 0.85*textSizeLabelsRel, 2, "", 42, 0.25);
  
  
  for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
    for (Int_t pid = 1; pid < nPID; pid++){
      if (!enableParticle[pid]) continue;
      if (debugOutput)std::cout << "effi E"  << std::endl;
      
      if (enableRecE){
        histoDummyEffiE->Draw();
        DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        legendEffiE->Clear();
        for(Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
//           if (!enablePlot[iEta]) continue;
          int iEtaPlot = iEta==maxEtaBinCaloDis[region] ? nEta : iEta;
          DrawGammaSetMarker(h_effi_rec_E[pid][iEta][iCl], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
          h_effi_rec_E[pid][iEta][iCl]->Draw("same,p");
          if ( h_effi_rec_E[pid][iEta][iCl]->GetMaximum() > 0){
            if (iEta == maxEtaBinCaloDis[region] )
              legendEffiE->AddEntry(h_effi_rec_E[pid][iEta][iCl],Form("%1.1f < #it{#eta} < %1.1f",partEta[minEtaBinCaloDis[region]],partEta[maxEtaBinCaloDis[region]]),"p");
            else 
              legendEffiE->AddEntry(h_effi_rec_E[pid][iEta][iCl],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
          }
        }
        legendEffiE->Draw();
        drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        if (nClusProcess != 1){
          drawLatexAdd(Form("%s in %s, %s clusters", partLabel[pid].Data(), caloPlot.Data(), nameClus[iCl].Data() ),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        } else {
           drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        cReso->Print(Form("%s/%s%s/Effi_E_%s.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), partName[pid].Data(), suffix.Data()));

        if (debugOutput)std::cout << "effi single entry E"  << std::endl;
        histoDummyEffiE->Draw();
        DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        for(Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
//           if (!enablePlot[iEta]) continue;
          int iEtaPlot = iEta==maxEtaBinCaloDis[region] ? nEta : iEta;
          DrawGammaSetMarker(h_effi_recSE_E[pid][iEta][iCl], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
          h_effi_recSE_E[pid][iEta][iCl]->Draw("same,p");
        }
        legendEffiE->Draw();
        drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        if (nClusProcess != 1){
          drawLatexAdd(Form("%s in %s, %s clusters", partLabel[pid].Data(), caloPlot.Data(), nameClus[iCl].Data() ),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        } else {
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), caloPlot.Data() ),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);          
        }
        if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        cReso->Print(Form("%s/%s%s/EffiSE_E_%s.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(),  partName[pid].Data(), suffix.Data()));
      }
      
      if (debugOutput)std::cout << "effi MCE"  << std::endl;
      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      legendEffiE->Clear();
      for(Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
//         if (!enablePlot[iEta]) continue;
        int iEtaPlot = iEta==maxEtaBinCaloDis[region] ? nEta : iEta;
        DrawGammaSetMarker(h_effi_rec_MCE[pid][iEta][iCl], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
        h_effi_rec_MCE[pid][iEta][iCl]->Draw("same,p");
        if ( h_effi_rec_MCE[pid][iEta][iCl]->GetMaximum() > 0){
          if (iEta == maxEtaBinCaloDis[region] )
            legendEffiE->AddEntry(h_effi_rec_MCE[pid][iEta][iCl],Form("%1.1f < #it{#eta} < %1.1f",partEta[minEtaBinCaloDis[region]],partEta[maxEtaBinCaloDis[region]]),"p");
          else 
            legendEffiE->AddEntry(h_effi_rec_MCE[pid][iEta][iCl],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        }
      }
      legendEffiE->Draw();
      drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      if (nClusProcess != 1){
        drawLatexAdd(Form("%s in %s, %s clusters", partLabel[pid].Data(), caloPlot.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      } else {
        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      }
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      cReso->Print(Form("%s/%s%s/Effi_MCE_%s.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), partName[pid].Data(),  suffix.Data()));

      if (debugOutput)std::cout << "effi single entry MCE"  << std::endl;
      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for(Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
//         if (!enablePlot[iEta]) continue;
        int iEtaPlot = iEta==maxEtaBinCaloDis[region] ? nEta : iEta;
        DrawGammaSetMarker(h_effi_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
        h_effi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");
      }
      legendEffiE->Draw();
      drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      if (nClusProcess != 1){
        drawLatexAdd(Form("%s in %s, %s clusters", partLabel[pid].Data(), caloPlot.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      } else {
        drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      }
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      cReso->Print(Form("%s/%s%s/EffiSE_MCE_%s.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), partName[pid].Data(),   suffix.Data()));

      if (enableTM){
        if (debugOutput)std::cout << "TM effi single entry"  << std::endl;
        histoDummyEffiTMMCE->Draw();
        DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        for(Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
//           if (!enablePlot[iEta]) continue;
          int iEtaPlot = iEta==maxEtaBinCaloDis[region] ? nEta : iEta;
          DrawGammaSetMarker(h_TMeffi_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
          h_TMeffi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");
        }
        legendEffiE->Draw();
        drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        if (nClusProcess != 1){
          drawLatexAdd(Form("%s in %s, %s clusters", partLabel[pid].Data(), caloPlot.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        } else {
           drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

        cReso->Print(Form("%s/%s%s/TMEffiSE_MCE_%s.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), partName[pid].Data(),   suffix.Data()));
  
        if (debugOutput)std::cout << "TM effi cluster-only single entry"  << std::endl;
        histoDummyEffiTMMCE->Draw();
        DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        for(Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
//           if (!enablePlot[iEta]) continue;
          int iEtaPlot = iEta==maxEtaBinCaloDis[region] ? nEta : iEta;
          DrawGammaSetMarker(h_TMeffiCls_recSE_MCE[pid][iEta][iCl], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
          h_TMeffiCls_recSE_MCE[pid][iEta][iCl]->Draw("same,p");
        }
        legendEffiE->Draw();
        drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        if (nClusProcess != 1){
          drawLatexAdd(Form("%s in %s, %s clusters", partLabel[pid].Data(), caloPlot.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        } else {
           drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

        cReso->Print(Form("%s/%s%s/TMEffiSEClus_MCE_%s.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), partName[pid].Data(),   suffix.Data()));
      }
    }
    for(Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
      Double_t etaMin = partEta[minEtaBinCaloDis[region]];
      Double_t etaMax = partEta[maxEtaBinCaloDis[region]];
      if (iEta < maxEtaBinCaloDis[region]){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }
      if (debugOutput)std::cout << "Effi PID MCE"  << std::endl;

      histoDummyEffiMCE->Draw();
      legendEffiPID->Clear();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t pid =1; pid < nPID; pid++){
        if (!enableParticle[pid]) continue;
        DrawGammaSetMarker(h_effi_rec_MCE[pid][iEta][iCl], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
        h_effi_rec_MCE[pid][iEta][iCl]->Draw("same,p");      
        legendEffiPID->AddEntry(h_effi_rec_MCE[pid][iEta][iCl],partLabel[pid].Data(),"p");
      }
      legendEffiPID->Draw();
      drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f < #it{#eta} < %1.1f, %s, %s clusters", etaMin, etaMax, caloPlot.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s%s/EffiPID_MCE_%d_%d.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      
      if (debugOutput)std::cout << "Effi PID single entry MCE"  << std::endl;
      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t pid =1; pid < nPID; pid++){
        if (!enableParticle[pid]) continue;
        DrawGammaSetMarker(h_effi_recSE_MCE[pid][iEta][iCl], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
        h_effi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
      }
      legendEffiPID->Draw();
      drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f < #it{#eta} < %1.1f, %s, %s clusters", etaMin, etaMax, caloPlot.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s%s/EffiPIDSE_MCE_%d_%d.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

      if (enableTM){
        if (debugOutput)std::cout << "TM Effi PID single entry MCE"  << std::endl;
        histoDummyEffiTMMCE->Draw();
        DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        for (Int_t pid =1; pid < nPID; pid++){
          if (!enableParticle[pid]) continue;
          DrawGammaSetMarker(h_TMeffi_recSE_MCE[pid][iEta][iCl], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
          h_TMeffi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
        }
        legendEffiPID->Draw();
        drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f < #it{#eta} < %1.1f, %s, %s clusters", etaMin, etaMax, caloPlot.Data(), nameClus[iCl].Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        cReso->Print(Form("%s/%s%s/TMEffiPIDSE_MCE_%d_%d.%s", outputDir.Data(), caloPlot.Data(), nameClus[iCl].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      }
    }
  }
  
  for (Int_t pid =1; pid < nPID; pid++){
    if (!enableParticle[pid]) continue;
    for (Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region];iEta++){
      Double_t etaMin = partEta[minEtaBinCaloDis[region]];
      Double_t etaMax = partEta[maxEtaBinCaloDis[region]];
      if (iEta < maxEtaBinCaloDis[region]){
        etaMin = partEta[iEta];
        etaMax = partEta[iEta+1];
      }
      if (debugOutput)std::cout << "Effi Clusterizer MCE"  << std::endl;
      histoDummyEffiMCE->Draw();
      legendEffiCl->Clear();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
        if (!enableClus[iCl]) continue;
        DrawGammaSetMarker(h_effi_rec_MCE[pid][iEta][iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
        h_effi_rec_MCE[pid][iEta][iCl]->Draw("same,p");      
        legendEffiCl->AddEntry(h_effi_rec_MCE[pid][iEta][iCl],nameClus[iCl].Data(),"p");
      }
      legendEffiCl->Draw();
      drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f < #it{#eta} < %1.1f, %s in %s", etaMin, etaMax, partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/EffiClusterizer_MCE_%s_%s_%d_%d.%s", outputDir.Data(), caloPlot.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      
      if (debugOutput)std::cout << "Effi Clusterizer single entry MCE"  << std::endl;
      histoDummyEffiMCE->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
      for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
        if (!enableClus[iCl]) continue;
        DrawGammaSetMarker(h_effi_recSE_MCE[pid][iEta][iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
        h_effi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
      }
      legendEffiCl->Draw();
      drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%1.1f < #it{#eta} < %1.1f, %s in %s", etaMin, etaMax, partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/EffiClusterizerSE_MCE_%s_%s_%d_%d.%s", outputDir.Data(), caloPlot.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

      if (enableTM){
        if (debugOutput)std::cout << "TM Effi Clusterizer single entry MCE"  << std::endl;
        histoDummyEffiTMMCE->Draw();
        DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
          if (!enableClus[iCl]) continue;
          DrawGammaSetMarker(h_TMeffi_recSE_MCE[pid][iEta][iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
          h_TMeffi_recSE_MCE[pid][iEta][iCl]->Draw("same,p");      
        }
        legendEffiCl->Draw();
        drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%1.1f < #it{#eta} < %1.1f, %s in %s", etaMin, etaMax, partLabel[pid].Data(), caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        cReso->Print(Form("%s/TMEffiClusterizerSE_MCE_%s_%s_%d_%d.%s", outputDir.Data(), caloPlot.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      }
    }
  }
  
  if (debugOutput)std::cout << "NTower"  << std::endl;
  TH2F* histoDummyNTowerMean   = new TH2F("histoDummyNTowerMean","histoDummyNTowerMean",1000,0, 200,1000,0.0, nMaxTowCl);
  SetStyleHistoTH2ForGraphs(histoDummyNTowerMean, "#it{E} (GeV)","#LT#it{N}_{tower}#GT", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNTowerMean->GetXaxis()->SetNoExponent();
  histoDummyNTowerMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNTowerMean->GetXaxis()->SetMoreLogLabels(kTRUE);

  TLegend* legendNTowerCl   = GetAndSetLegend2(0.12,  0.94-(nActiceCl/2*0.85*textSizeLabelsRel), 0.4,  0.94, 0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  histoDummyNTowerMean->Draw();
  for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
    if (!enableClus[iCl]) continue;
    DrawGammaSetMarker(h_cluster_NTowerMean_E[iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
    h_cluster_NTowerMean_E[iCl]->Draw("same,p");      
    legendNTowerCl->AddEntry(h_cluster_NTowerMean_E[iCl],nameClus[iCl].Data(),"p");
  }
  legendNTowerCl->Draw();
  drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("%s", caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/NTowerInCluster_%s_Mean_E.%s", outputDir.Data(), caloPlot.Data(),  suffix.Data()));

  
  
  if (debugOutput)std::cout << "NCluster"  << std::endl;
  TH2F* histoDummyNClPerParMean   = new TH2F("histoDummyNClPerParMean","histoDummyNClPerParMean",1000,0, 200,1000,0.0, nMaxClPart);
  SetStyleHistoTH2ForGraphs(histoDummyNClPerParMean, "#it{E} (GeV)","#LT#it{N}_{cl}/particle #GT", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNClPerParMean->GetXaxis()->SetNoExponent();
  histoDummyNClPerParMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNClPerParMean->GetXaxis()->SetMoreLogLabels(kTRUE);

  if (enableRecE){
    histoDummyNClPerParMean->Draw();
    DrawGammaLines(0.2, 200, 1., 1., 2, kGray+2, 7);
    for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
      if (!enableClus[iCl]) continue;
      DrawGammaSetMarker(h_cluster_NClMean_E[iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
      h_cluster_NClMean_E[iCl]->Draw("same,p");      
    }
    legendNTowerCl->Draw();
    drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s", caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/NClusterPerParticle_%s_Mean_E.%s", outputDir.Data(), caloPlot.Data(),  suffix.Data()));
  }
    
  if (debugOutput)std::cout << "NCluster MCE"  << std::endl;
  histoDummyNClPerParMean->SetXTitle("#it{E}_{MC} (GeV)");
  histoDummyNClPerParMean->Draw();
  DrawGammaLines(0.2, 200, 1., 1., 2, kGray+2, 7);
  for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
    if (!enableClus[iCl]) continue;
    DrawGammaSetMarker(h_cluster_NClMean_MCE[iCl], markerStyleClus[iCl], markerSizeClus[iCl], colorClus[0][iCl], colorClus[0][iCl]);
    h_cluster_NClMean_MCE[iCl]->Draw("same,p");      
  }
  legendNTowerCl->Draw();
  drawLatexAdd(labelEnergy,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("%s", caloPlot.Data()),0.95,0.91-nLinesCol*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/NClusterPerParticle_%s_Mean_MCE.%s", outputDir.Data(), caloPlot.Data(),  suffix.Data()));
  
//    
  TFile* outputFile  = new TFile(inputFileNameCluster.Data(),"UPDATE");
  TDirectoryFile* directoryClEffi = (TDirectoryFile*)outputFile->Get(Form("ClusterEfficiencyProjections_%s", calo.Data()));
  if (!directoryClEffi){
    outputFile->mkdir(Form("ClusterEfficiencyProjections_%s", calo.Data()));
    directoryClEffi = (TDirectoryFile*)outputFile->Get(Form("ClusterEfficiencyProjections_%s", calo.Data()));
  }
  directoryClEffi->cd();
  for (Int_t iCl = 0; iCl < nClusProcess; iCl++){
    for (Int_t iEta=minEtaBinCaloDis[region]; iEta<maxEtaBinCaloDis[region]+1;iEta++){
      for (Int_t pid = 1; pid < 6; pid++){
        if (!enableParticle[pid]) continue;
        if (h_effi_rec_E[pid][iEta][iCl]) h_effi_rec_E[pid][iEta][iCl]->Write(Form("effi%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        if (h_effi_rec_MCE[pid][iEta][iCl]) h_effi_rec_MCE[pid][iEta][iCl]->Write(Form("effi%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        if (h_effi_recSE_E[pid][iEta][iCl]) h_effi_recSE_E[pid][iEta][iCl]->Write(Form("effiSE%s_E_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        if (h_effi_recSE_MCE[pid][iEta][iCl]) h_effi_recSE_MCE[pid][iEta][iCl]->Write(Form("effiSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        if (h_TMeffi_recSE_MCE[pid][iEta][iCl]) h_TMeffi_recSE_MCE[pid][iEta][iCl]->Write(Form("TMeffiSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
        if (h_TMeffiCls_recSE_MCE[pid][iEta][iCl]) h_TMeffiCls_recSE_MCE[pid][iEta][iCl]->Write(Form("TMeffiClsSE%s_MCE_%d_%s",partName[pid].Data(), iEta, nameClus[iCl].Data()),TObject::kOverwrite);
      }
    }
    h_cluster_NTowerMean_E[iCl]->Write(Form("h_CS_NTowerMean_%s_E", nameClus[iCl].Data()),TObject::kOverwrite);
    if(h_cluster_NClMean_E[iCl]) h_cluster_NClMean_E[iCl]->Write(Form("h_CS_NClMean_%s_E", nameClus[iCl].Data()),TObject::kOverwrite);
    h_cluster_NClMean_MCE[iCl]->Write(Form("h_CS_NClMean_%s_MCE", nameClus[iCl].Data()),TObject::kOverwrite);
  }
  outputFile->Write();
  outputFile->Close();
}
