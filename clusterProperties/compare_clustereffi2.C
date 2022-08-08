#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void MaskHighEPoint (TH1D* histo, Double_t maxE = 1000, bool mask = false, Double_t value = 0.){
  if (mask && histo){
    for (Int_t iE = 1; iE < histo->GetNbinsX()+1; iE++){
      if (histo->GetBinCenter(iE) > maxE){
        histo->SetBinContent(iE, value); 
        histo->SetBinError(iE, 0); 
      }
      if (histo->GetBinCenter(iE) < 0.2){
        histo->SetBinContent(iE, value); 
        histo->SetBinError(iE, 0); 
      }
    }
  }
}


void compare_clustereffi2(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString suffix            = "pdf",
                            bool paperPlot            = true
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();
  TString detLabel = "";

  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;
  Double_t maxBetaSigma       = 0.05;

  TString outputDir                 = Form("plots/%s/Compare%s",dateForOutput.Data(), addName.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TString labelEnergy   = "#it{#bf{ECCE}} simulation";
  TString collisionSystem = GetCollisionEnergy(addName);
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addName.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }
    
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  Color_t colorSet[8]         = {kBlack, kRed+1, kGreen+2, kGray+1, kCyan+1, kAzure+2, kBlack };
  Style_t markerStyleSet[8]   = {20, 47, 33,  21, 30, 42, 46, 20};
  Size_t markerSizeSet[8]     = {1.8, 2.0, 2.4, 1.5, 2.4, 1.8, 1.5, 1.5 };
  
  const Int_t maxNSets        = 10;
  TString outTrackCuts[3]     = {"", "LI2", "LI3"};
  TString labelTrackCuts[3]   = {"", "#geq 2 tracker hits", "#geq 3 tracker hits"};

  TString outBetaCuts[2]      = {"", "LI3"};
  TString labelBetaCuts[2]    = {"", "#geq 3 tracker hits"};
  
  
  TH1D* h_effi_rec_E[maxNSets][nPID][nEta+1]          = {{{NULL}}};
  TH1D* h_effi_rec_MCE[maxNSets][nPID][nEta+1]        = {{{NULL}}};
  TH1D* h_effi_recSE_E[maxNSets][nPID][nEta+1]        = {{{NULL}}};
  TH1D* h_effi_recSE_MCE[maxNSets][nPID][nEta+1]      = {{{NULL}}};
  TH1D* h_TMeffi_recSE_MCE[maxNSets][nPID][nEta+1]    = {{{NULL}}};
  TH1D* h_TMeffiCls_recSE_MCE[maxNSets][nPID][nEta+1] = {{{NULL}}};
  TH1D* hSP_effi_rec_E[maxNSets][nPID][3]             = {{{NULL}}};
  TH1D* hSP_effi_rec_MCE[maxNSets][nPID][3]           = {{{NULL}}};
  TH1D* hSP_effi_recSE_E[maxNSets][nPID][3]           = {{{NULL}}};
  TH1D* hSP_effi_recSE_MCE[maxNSets][nPID][3]         = {{{NULL}}};
  TH1D* hSP_TMeffi_recSE_MCE[maxNSets][nPID][3]       = {{{NULL}}};
  TH1D* hSP_TMeffiCls_recSE_MCE[maxNSets][nPID][3]    = {{{NULL}}};
  TH1D* h_cluster_NTowerMean_E[maxNSets]              = {NULL};
  TH1D* h_cluster_NClMean_E[maxNSets]                 = {NULL};
  TH1D* h_cluster_NClMean_MCE[maxNSets]               = {NULL};
  TFile* inputFiles[maxNSets]                         = {NULL};
  TString inputFilesNames[maxNSets];
  TString labels[maxNSets];
  TString calo[maxNSets];
  int region[maxNSets];
  int caloID[maxNSets];
  Double_t nomEtaRange[maxNSets][2]     = {{0.}} ;
  Double_t maxEEcal = 30;
  Double_t maxHEcal = 50;
  
  bool currPIDavail[nPID] = {false};

  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> calo[nSets] >> region[nSets]>> labels[nSets] >> caloID[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< calo[nSets].Data()<< "\t"<< labels[nSets].Data() << "\t" << caloID[nSets] <<std::endl;
      labels[nSets].ReplaceAll("_"," ");
      colorSet[nSets] = getCaloColor(calo[nSets]);
      nomEtaRange[nSets][0] = nominalEtaRegion[region[nSets]][0];
      nomEtaRange[nSets][1] = nominalEtaRegion[region[nSets]][1];
      if (caloID[nSets] == 4) // iHCAL
        nomEtaRange[nSets][0] = -1.3;
      else if (caloID[nSets] == 5) // oHCAL
        nomEtaRange[nSets][0] = -1.1;
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;


  for (Int_t iSet = 0; iSet < nSets; iSet++){
    colorSet[iSet]        = getCaloColor(labels[iSet],false);
    markerStyleSet[iSet]  = getCaloMarker(labels[iSet]);
    inputFiles[iSet]  = new TFile(inputFilesNames[iSet].Data());
    for (Int_t iEta=minEtaBinCaloDis[region[iSet]]; iEta<maxEtaBinCaloDis[region[iSet]]+1;iEta++){
      for (Int_t pid = 1; pid < nPID; pid++){
//             cout << Form("ClusterEfficiencyProjections_%s/effi%s_E_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA") << endl;
        h_effi_rec_E[iSet][pid][iEta]      = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effi%s_E_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        h_effi_rec_MCE[iSet][pid][iEta]    = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effi%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        h_effi_recSE_E[iSet][pid][iEta]    = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiSE%s_E_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        h_effi_recSE_MCE[iSet][pid][iEta]  = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        if(h_effi_recSE_MCE[iSet][pid][iEta]) currPIDavail[pid] = true;
        else cout << Form("ClusterEfficiencyProjections_%s/effiSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA") << " not found" << endl;  
        h_TMeffi_recSE_MCE[iSet][pid][iEta]  = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/TMeffiSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        h_TMeffiCls_recSE_MCE[iSet][pid][iEta]  = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/TMeffiClsSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        if (currPIDavail[pid]){
          if (iSet < 3){
            MaskHighEPoint(h_effi_rec_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
            MaskHighEPoint(h_effi_recSE_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
            MaskHighEPoint(h_TMeffi_recSE_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
            MaskHighEPoint(h_TMeffiCls_recSE_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
          } else {
            MaskHighEPoint(h_effi_rec_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
            MaskHighEPoint(h_effi_recSE_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
            MaskHighEPoint(h_TMeffi_recSE_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
            MaskHighEPoint(h_TMeffiCls_recSE_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
          }          
        }
      }
    }
    for (Int_t iEta=0; iEta<3;iEta++){
      for (Int_t pid = 1; pid < nPID; pid++){
        hSP_effi_rec_E[iSet][pid][iEta]           = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiPaper%s_E_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        hSP_effi_rec_MCE[iSet][pid][iEta]         = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiPaper%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        hSP_effi_recSE_E[iSet][pid][iEta]         = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiPaperSE%s_E_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        hSP_effi_recSE_MCE[iSet][pid][iEta]       = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiPaperSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        if(hSP_effi_recSE_MCE[iSet][pid][iEta]) currPIDavail[pid] = true;
        else cout << Form("ClusterEfficiencyProjections_%s/effiPaperSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA") << " not found" << endl;
        hSP_TMeffi_recSE_MCE[iSet][pid][iEta]     = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/TMeffiPaperSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        hSP_TMeffiCls_recSE_MCE[iSet][pid][iEta]  = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/TMeffiPaperClsSE%s_MCE_%d_%s", calo[iSet].Data(), partName[pid].Data(), iEta, "MA"));
        if (currPIDavail[pid]){
          if (iSet < 3){
            MaskHighEPoint(hSP_effi_rec_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
            MaskHighEPoint(hSP_effi_recSE_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
            MaskHighEPoint(hSP_TMeffi_recSE_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
            MaskHighEPoint(hSP_TMeffiCls_recSE_MCE[iSet][pid][iEta], maxEEcal, paperPlot);
          } else {
            MaskHighEPoint(hSP_effi_rec_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
            MaskHighEPoint(hSP_effi_recSE_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
            MaskHighEPoint(hSP_TMeffi_recSE_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
            MaskHighEPoint(hSP_TMeffiCls_recSE_MCE[iSet][pid][iEta], maxHEcal, paperPlot);
          }          
        }
      }
    }
    h_cluster_NTowerMean_E[iSet] = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NTowerMean_%s_E", calo[iSet].Data(), "MA"));
    h_cluster_NClMean_E[iSet]    = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NClMean_%s_E", calo[iSet].Data(), "MA"));
    h_cluster_NClMean_MCE[iSet]  = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NClMean_%s_MCE", calo[iSet].Data(), "MA"));
    if (iSet < 3){
      MaskHighEPoint(h_cluster_NTowerMean_E[iSet], maxEEcal, paperPlot);
      MaskHighEPoint(h_cluster_NClMean_E[iSet], maxEEcal, paperPlot);
      MaskHighEPoint(h_cluster_NClMean_MCE[iSet], maxEEcal, paperPlot);
    } else {
      MaskHighEPoint(h_cluster_NTowerMean_E[iSet], maxHEcal, paperPlot);
      MaskHighEPoint(h_cluster_NClMean_E[iSet], maxHEcal, paperPlot);
      MaskHighEPoint(h_cluster_NClMean_MCE[iSet], maxHEcal, paperPlot);    
    }
    
  }

  TCanvas* cReso = new TCanvas("cReso","",0,0,900,800);
  DrawGammaCanvasSettings( cReso, 0.085, 0.025, 0.02, 0.105);
  cReso->SetLogx();

  TH2F* histoDummyEffiE   = new TH2F("histoDummyEffiE","histoDummyEffiE",1000,0.15, 70,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiE, "#it{E} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.87);
  histoDummyEffiE->GetXaxis()->SetNoExponent();
  histoDummyEffiE->GetXaxis()->SetRangeUser(0.15,50);
  histoDummyEffiE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiE->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 70,1000,0.0, 1.55);
  SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.57);
  histoDummyEffiMCE->GetXaxis()->SetNoExponent();
  histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
  histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,50);
  histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);

  TH2F* histoDummyNTowerMean   = new TH2F("histoDummyNTowerMean","histoDummyNTowerMean",1000,0.21, 70,1000,1, 999);
  SetStyleHistoTH2ForGraphs(histoDummyNTowerMean, "#it{E} (GeV)","#LT#it{N}_{tower}#GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNTowerMean->GetXaxis()->SetNoExponent();
  histoDummyNTowerMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNTowerMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyNClPerParMean   = new TH2F("histoDummyNClPerParMean","histoDummyNClPerParMean",1000,0.21, 70,1000,0.0, 4.9);
  SetStyleHistoTH2ForGraphs(histoDummyNClPerParMean, "#it{E} (GeV)","#LT#it{N}_{cl}/particle #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNClPerParMean->GetXaxis()->SetNoExponent();
  histoDummyNClPerParMean->GetYaxis()->SetNdivisions(505,kTRUE);
  histoDummyNClPerParMean->GetXaxis()->SetMoreLogLabels(kTRUE);


  TLegend* legendPtResMLeft  = GetAndSetLegend2(0.12, 0.82-(nSets*textSizeLabelsRel), 0.30, 0.82,textSizeLabelsPixel, 1, "", 43, 0.15);
  TLegend* legendPtResM  = GetAndSetLegend2(0.74, 0.935-(nSets*textSizeLabelsRel), 0.92, 0.935,textSizeLabelsPixel, 1, "", 43, 0.25);
  legendPtResM->SetMargin(0.3);
  legendPtResMLeft->SetMargin(0.3);
  for (Int_t pid = 1; pid < nPID; pid++){
    if(!currPIDavail[pid]) continue;
      legendPtResM->Clear();
    histoDummyNClPerParMean->Draw();
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_cluster_NClMean_E[iSet], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_cluster_NClMean_E[iSet]->Draw("same,hist,p");
        legendPtResM->AddEntry(h_cluster_NClMean_E[iSet],labels[iSet].Data(),"p");
        legendPtResMLeft->AddEntry(h_cluster_NClMean_E[iSet],labels[iSet].Data(),"p");
      }
      legendPtResM->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single particles"),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cReso->Print(Form("%s/%s_NClusterPerParticle_Mean_MCE.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    
  cReso->SetLogy(true);    
    histoDummyNTowerMean->Draw();
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_cluster_NTowerMean_E[iSet], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_cluster_NTowerMean_E[iSet]->Draw("same,hist,p");
      }
      legendPtResMLeft->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single particles"),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print(Form("%s/%s_NTowerInCluster_Mean_E.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
  cReso->SetLogy(false);
    
    histoDummyEffiMCE->Draw();
    DrawGammaLines(0.15, 70, 1,1,2,kGray+1, 2);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]]->Draw("same,p");
      }
      legendPtResM->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single %s", partLabel[pid].Data()),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print(Form("%s/%s_ClusterTMEfficiencySE_MCE.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    
    histoDummyEffiMCE->Draw();
    DrawGammaLines(0.15, 70, 1,1,2,kGray+1, 2);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_TMeffiCls_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_TMeffiCls_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]]->Draw("same,p");
        // legendPtResM->AddEntry(h_TMeffiCls_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]],labels[iSet].Data(),"p");
      }
      legendPtResM->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single %s", partLabel[pid].Data()),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{clus}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print(Form("%s/%s_ClusterTMEfficiencyClusterSE_MCE.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
  }
  
  //*********************************************************************************************************
  //*********************************************************************************************************
  //=============================== PAPER PLOTTING ECCE NIM CALO ============================================
  //*********************************************************************************************************
  //*********************************************************************************************************
  
  if(nSets==5 && paperPlot){

    //*********************************************************************************************************
    //***** Cluster property paper plots, top mean cluster per particle, bottom mean cells per cluster ********
    //*********************************************************************************************************
    Double_t arrayBound2PanY_X[2];
    Double_t arrayBound2PanY_Y[3];
    Double_t relMargin2PanY_X[3];
    Double_t relMargin2PanY_Y[3];
    textSizeLabelsPixel             = 50;
    ReturnCorrectValuesForCanvasScaling(1350,1450, 1, 2,0.07, 0.006, 0.005,0.075,arrayBound2PanY_X,arrayBound2PanY_Y,relMargin2PanY_X,relMargin2PanY_Y);

    TCanvas* canvasClusProp     = new TCanvas("canvasClusProp","",0,0,1350,1450);  // gives the page size
    DrawGammaCanvasSettings( canvasClusProp,  0.13, 0.02, 0.03, 0.06);

    TPad* padClusPropTop               = new TPad("padClusPropTop", "", arrayBound2PanY_X[0], arrayBound2PanY_Y[1], arrayBound2PanY_X[1], arrayBound2PanY_Y[0],-1, -1, -2);
    DrawGammaPadSettings( padClusPropTop, relMargin2PanY_X[0], relMargin2PanY_X[2], relMargin2PanY_Y[0], relMargin2PanY_Y[1]);
    padClusPropTop->Draw();

    TPad* padClusPropBottom                = new TPad("padClusPropBottom", "", arrayBound2PanY_X[0], arrayBound2PanY_Y[2], arrayBound2PanY_X[1], arrayBound2PanY_Y[1],-1, -1, -2);
    DrawGammaPadSettings( padClusPropBottom, relMargin2PanY_X[0], relMargin2PanY_X[2], relMargin2PanY_Y[1], relMargin2PanY_Y[2]);
    padClusPropBottom->Draw();

    padClusPropTop->cd();
    padClusPropTop->SetLogx();

      Double_t textsizeLabelsNClMean    = 0;
      if (padClusPropTop->XtoPixel(padClusPropTop->GetX2()) < padClusPropTop->YtoPixel(padClusPropTop->GetY1())){
          textsizeLabelsNClMean         = (Double_t)textSizeLabelsPixel/padClusPropTop->XtoPixel(padClusPropTop->GetX2()) ;
      } else {
          textsizeLabelsNClMean         = (Double_t)textSizeLabelsPixel/padClusPropTop->YtoPixel(padClusPropTop->GetY1());
      }

      TH2F * histo2DNClMean    = new TH2F("histo2DNClMean","histo2DNClMean",1000,0.15, 70,1000,0.8, 2.8);
      SetStyleHistoTH2ForGraphs(histo2DNClMean, "#it{E} (GeV)","#LT#it{N}_{cl}/particle#GT", 0.85*textsizeLabelsNClMean, textsizeLabelsNClMean,
                                0.85*textsizeLabelsNClMean, textsizeLabelsNClMean, 0.8,0.35,512,502);//#it{p}_{T} (GeV/#it{c})
      histo2DNClMean->GetYaxis()->SetMoreLogLabels(kTRUE);
      histo2DNClMean->GetYaxis()->SetNoExponent(kTRUE);
      histo2DNClMean->GetXaxis()->SetTickLength(0.05);
      histo2DNClMean->GetYaxis()->SetTickLength(0.026);
      histo2DNClMean->GetYaxis()->SetLabelOffset(0.005);
      histo2DNClMean->SetMaximum(20);
      histo2DNClMean->DrawCopy();

      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_cluster_NClMean_E[iSet], markerStyleSet[iSet], 1.5*markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_cluster_NClMean_E[iSet]->Draw("same,hist,p");
      }  
      drawLatexAdd("#it{#bf{ECCE}} simulation",0.10,0.86,textsizeLabelsNClMean,kFALSE,kFALSE,false);
      drawLatexAdd(Form("single particles"),0.10,0.78,textsizeLabelsNClMean,kFALSE,kFALSE,false);
      drawLatexAdd(Form("a)"),0.95,0.08,textsizeLabelsNClMean,kFALSE,kFALSE,true);
    padClusPropBottom->SetLogx();
    padClusPropBottom->SetLogy();

    padClusPropBottom->cd();
      Double_t textsizeLabelsNtowMean         = 0;
      if (padClusPropBottom->XtoPixel(padClusPropBottom->GetX2()) <padClusPropBottom->YtoPixel(padClusPropBottom->GetY1()) ){
          textsizeLabelsNtowMean              = (Double_t)textSizeLabelsPixel/padClusPropBottom->XtoPixel(padClusPropBottom->GetX2()) ;
      } else {
          textsizeLabelsNtowMean              = (Double_t)textSizeLabelsPixel/padClusPropBottom->YtoPixel(padClusPropBottom->GetY1());
      }

      TH2F * histo2DNTowerMean            = new TH2F("histo2DNTowerMean","histo2DNTowerMean",1000,0.15, 70,1000,1, 999);//, 100.1, 160.9);//125.1, 155.9);
      SetStyleHistoTH2ForGraphs(histo2DNTowerMean, "#it{E} (GeV)","#LT#it{N}_{tower}#GT", 0.85*textsizeLabelsNtowMean, textsizeLabelsNtowMean, 0.85*textsizeLabelsNtowMean,
                                textsizeLabelsNtowMean, 0.9, 0.35,512,505);
      histo2DNTowerMean->GetXaxis()->SetMoreLogLabels(kTRUE);
      histo2DNTowerMean->GetYaxis()->SetNoExponent(kTRUE);
      histo2DNTowerMean->GetXaxis()->SetTickLength(0.05);
      histo2DNTowerMean->GetXaxis()->SetNoExponent();
      histo2DNTowerMean->SetMaximum(20);
      histo2DNTowerMean->DrawCopy();
      TLegend* legendPtResMLeft2  = GetAndSetLegend2(0.10, 0.95-(nSets*textsizeLabelsNtowMean), 0.28, 0.95,textSizeLabelsPixel, 1, "", 43, 0.18);

      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_cluster_NTowerMean_E[iSet], markerStyleSet[iSet], 1.5*markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
        h_cluster_NTowerMean_E[iSet]->Draw("same,hist,p");
        legendPtResMLeft2->AddEntry(h_cluster_NTowerMean_E[iSet],labels[iSet].Data(),"p");
      }
      legendPtResMLeft2->Draw();
      drawLatexAdd(Form("b)"),0.95,0.19,textsizeLabelsNtowMean,kFALSE,kFALSE,true);

    canvasClusProp->Update();
    canvasClusProp->Print(Form("%s/NTowAndMeanClus_Paper.%s",outputDir.Data(),suffix.Data()));

    
    //*********************************************************************************************************
    //************************************ TM effi summary plot  **********************************************
    //*********************************************************************************************************        
    padClusPropTop->cd();
    padClusPropTop->SetLogx();

      Double_t textsizeLabelsEffiTMTop    = 0;
      if (padClusPropTop->XtoPixel(padClusPropTop->GetX2()) < padClusPropTop->YtoPixel(padClusPropTop->GetY1())){
          textsizeLabelsEffiTMTop         = (Double_t)textSizeLabelsPixel/padClusPropTop->XtoPixel(padClusPropTop->GetX2()) ;
      } else {
          textsizeLabelsEffiTMTop         = (Double_t)textSizeLabelsPixel/padClusPropTop->YtoPixel(padClusPropTop->GetY1());
      }

      TH2F * histo2DAllEffiTM2    = new TH2F("histo2DAllEffiTM2","histo2DAllEffiTM2",1000,0.15, 70,1000,0.001, 1.35);
      SetStyleHistoTH2ForGraphs(histo2DAllEffiTM2, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textsizeLabelsEffiTMTop, textsizeLabelsEffiTMTop,
                                0.85*textsizeLabelsEffiTMTop, textsizeLabelsEffiTMTop, 0.8,0.45,512,505);//#it{p}_{T} (GeV/#it{c})
      histo2DAllEffiTM2->GetYaxis()->SetMoreLogLabels(kTRUE);
      histo2DAllEffiTM2->GetYaxis()->SetNoExponent(kTRUE);
      histo2DAllEffiTM2->GetXaxis()->SetTickLength(0.05);
      histo2DAllEffiTM2->GetYaxis()->SetTickLength(0.026);
      histo2DAllEffiTM2->SetMaximum(20);
      histo2DAllEffiTM2->DrawCopy();
      DrawGammaLines(0.15, 70, 1,1,2,kGray+1, 2);
        for (Int_t pid = 1; pid < nPID; pid++){
        if(!currPIDavail[pid]) continue;
        for(Int_t iSet=0; iSet<nSets;iSet++){
            DrawGammaSetMarker(h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], getCaloMarker(labels[iSet]), 1.5*markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        }
      }
      h_TMeffi_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]]->Draw("same,p");

      drawLatexAdd("#it{#bf{ECCE}} simulation",0.10,0.86,textsizeLabelsEffiTMTop,kFALSE,kFALSE,false);
      drawLatexAdd(Form("single particles"),0.10,0.78,textsizeLabelsEffiTMTop,kFALSE,kFALSE,false);
      drawLatexAdd(Form("a)"),0.97,0.08,textsizeLabelsEffiTMTop,kFALSE,kFALSE,true);
      histo2DAllEffiTM2->Draw("same,axis");
    padClusPropBottom->SetLogx();
    padClusPropBottom->SetLogy(false);

      padClusPropBottom->cd();
      Double_t textsizeLabelsEffiTMBottom         = 0;
      if (padClusPropBottom->XtoPixel(padClusPropBottom->GetX2()) <padClusPropBottom->YtoPixel(padClusPropBottom->GetY1()) ){
          textsizeLabelsEffiTMBottom              = (Double_t)textSizeLabelsPixel/padClusPropBottom->XtoPixel(padClusPropBottom->GetX2()) ;
      } else {
          textsizeLabelsEffiTMBottom              = (Double_t)textSizeLabelsPixel/padClusPropBottom->YtoPixel(padClusPropBottom->GetY1());
      }

      TH2F * histo2DAllEffiTM3            = new TH2F("histo2DAllEffiTM3","histo2DAllEffiTM3",1000,0.15, 70,1000,0.0, 1.35);//, 100.1, 160.9);//125.1, 155.9);
      SetStyleHistoTH2ForGraphs(histo2DAllEffiTM3, "#it{E}^{MC} (GeV)","#kappa_{TM}", 0.85*textsizeLabelsEffiTMBottom, textsizeLabelsEffiTMBottom, 0.85*textsizeLabelsEffiTMBottom,
                                textsizeLabelsEffiTMBottom, 0.9, 0.45,512,505);
      histo2DAllEffiTM3->GetXaxis()->SetMoreLogLabels(kTRUE);
      histo2DAllEffiTM3->GetYaxis()->SetNoExponent(kTRUE);
      histo2DAllEffiTM3->GetXaxis()->SetTickLength(0.05);
      histo2DAllEffiTM3->GetXaxis()->SetNoExponent();
      histo2DAllEffiTM3->SetMaximum(20);
      histo2DAllEffiTM3->DrawCopy();
      legendPtResMLeft2  = GetAndSetLegend2(0.10, 0.94-(2*textsizeLabelsEffiTMBottom), 0.78, 0.94,textSizeLabelsPixel, 3, "", 43, 0.15);
      for (Int_t pid = 1; pid < nPID; pid++){
        if(!currPIDavail[pid]) continue;
        for(Int_t iSet=0; iSet<nSets;iSet++){
            DrawGammaSetMarker(h_TMeffiCls_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], getCaloMarker(labels[iSet]), 1.5*markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        }
      }
      DrawGammaLines(0.15, 70, 1,1,2,kGray+1, 2);
      h_TMeffiCls_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]]->Draw("same,p");

      legendPtResMLeft2->AddEntry(h_TMeffiCls_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]],Form("%s (%s)",labels[0].Data(),partLabel[1].Data()),"p");
      legendPtResMLeft2->AddEntry(h_TMeffiCls_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]],Form("%s (%s)",labels[1].Data(),partLabel[1].Data()),"p");
      legendPtResMLeft2->AddEntry(h_TMeffiCls_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]],Form("%s (%s)",labels[2].Data(),partLabel[1].Data()),"p");
      legendPtResMLeft2->AddEntry(h_TMeffiCls_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]],Form("%s (%s)",labels[3].Data(),partLabel[3].Data()),"p");
      legendPtResMLeft2->AddEntry(h_TMeffiCls_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]],Form("%s (%s)",labels[4].Data(),partLabel[3].Data()),"p");

      legendPtResMLeft2->Draw();
      drawLatexAdd(Form("b)"),0.97,0.19,textsizeLabelsEffiTMBottom,kFALSE,kFALSE,true);
      histo2DAllEffiTM3->Draw("same,axis");

    canvasClusProp->Update();
    canvasClusProp->Print(Form("%s/ClusterTMEfficiency_Paper.%s",outputDir.Data(),suffix.Data()));


    TCanvas* canvasTMEffKappa     = new TCanvas("canvasTMEffKappa","",0,0,1350,950);  // gives the page size
    DrawGammaCanvasSettings( canvasTMEffKappa,  0.07, 0.006, 0.005,0.105);
    canvasTMEffKappa->SetLogx();
      Double_t textsizeLabelsKappa    = 0;
      if (canvasTMEffKappa->XtoPixel(canvasTMEffKappa->GetX2()) < canvasTMEffKappa->YtoPixel(canvasTMEffKappa->GetY1())){
          textsizeLabelsKappa         = (Double_t)textSizeLabelsPixel/canvasTMEffKappa->XtoPixel(canvasTMEffKappa->GetX2()) ;
      } else {
          textsizeLabelsKappa         = (Double_t)textSizeLabelsPixel/canvasTMEffKappa->YtoPixel(canvasTMEffKappa->GetY1());
      }

      TH2F * histo2DAllEffiTM4            = new TH2F("histo2DAllEffiTM4","histo2DAllEffiTM4",1000,0.15, 70,1000,0.0, 1.35);//, 100.1, 160.9);//125.1, 155.9);
      SetStyleHistoTH2ForGraphs(histo2DAllEffiTM4, "#it{E}^{MC} (GeV)","#kappa_{TM}", 0.85*textsizeLabelsKappa, textsizeLabelsKappa, 0.85*textsizeLabelsKappa,
                                textsizeLabelsKappa, 0.9, 0.6,512,505);
      histo2DAllEffiTM4->GetXaxis()->SetMoreLogLabels(kTRUE);
      histo2DAllEffiTM4->GetYaxis()->SetNoExponent(kTRUE);
      histo2DAllEffiTM4->GetXaxis()->SetNoExponent();
      histo2DAllEffiTM4->SetMaximum(20);
      histo2DAllEffiTM4->DrawCopy();
      TLegend* legendPtResMLeft3  = GetAndSetLegend2(0.10, 0.94-(2*textsizeLabelsKappa), 0.78, 0.94,textSizeLabelsPixel, 3, "", 43, 0.15);
      for (Int_t pid = 1; pid < nPID; pid++){
        if(!currPIDavail[pid]) continue;
        for(Int_t iSet=0; iSet<nSets;iSet++){
            DrawGammaSetMarker(h_TMeffiCls_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], getCaloMarker(labels[iSet]), 1.5*markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        }
      }
      DrawGammaLines(0.15, 70, 1,1,2,kGray+1, 2);
      h_TMeffiCls_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]]->Draw("same,p");
      h_TMeffiCls_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]]->Draw("same,p");

      legendPtResMLeft3->AddEntry(h_TMeffiCls_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]],Form("%s (%s)",labels[0].Data(),partLabel[1].Data()),"p");
      legendPtResMLeft3->AddEntry(h_TMeffiCls_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]],Form("%s (%s)",labels[1].Data(),partLabel[1].Data()),"p");
      legendPtResMLeft3->AddEntry(h_TMeffiCls_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]],Form("%s (%s)",labels[2].Data(),partLabel[1].Data()),"p");
      legendPtResMLeft3->AddEntry(h_TMeffiCls_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]],Form("%s (%s)",labels[3].Data(),partLabel[3].Data()),"p");
      legendPtResMLeft3->AddEntry(h_TMeffiCls_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]],Form("%s (%s)",labels[4].Data(),partLabel[3].Data()),"p");

      legendPtResMLeft3->Draw();
      drawLatexAdd(labelEnergy,0.95,0.18+2*textsizeLabelsKappa*1.1,textsizeLabelsKappa,kFALSE,kFALSE,true);
      drawLatexAdd(Form("single particles"),0.95,0.18+textsizeLabelsKappa*1.1,textsizeLabelsKappa,kFALSE,kFALSE,true);
      drawLatexAdd(Form("#kappa_{TM} = N_{clus}^{matched} / N_{clus}"),0.95,0.17,textsizeLabelsKappa,kFALSE,kFALSE,true);
      histo2DAllEffiTM4->Draw("same,axis");

    canvasTMEffKappa->Update();
    canvasTMEffKappa->Print(Form("%s/ClusterTMEfficiency_KappaOnly_Paper.%s",outputDir.Data(),suffix.Data()));
        
    //************************************************************************************************
    //************ kappa: TM effi plot ECals w/ respect to found clusters all eta regions ************
    //************************************************************************************************  
    cReso->cd();
    histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.55);
      histoDummyEffiMCE->Draw();
      legendPtResM  = GetAndSetLegend2(0.68, 0.935-(nSets*textSizeLabelsRel), 0.86, 0.935,textSizeLabelsPixel, 1, "", 43, 0.25);
      legendPtResM->SetMargin(0.3);
      DrawGammaLines(0.15, 70, 1,1,2,kGray+1, 2);
      h_TMeffi_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]]->Draw("same,p");
      h_TMeffi_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]]->Draw("same,p");
      legendPtResM->AddEntry(h_TMeffi_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]],Form("%s (%s)",labels[0].Data(),partLabel[1].Data()),"p");
      legendPtResM->AddEntry(h_TMeffi_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]],Form("%s (%s)",labels[1].Data(),partLabel[1].Data()),"p");
      legendPtResM->AddEntry(h_TMeffi_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]],Form("%s (%s)",labels[2].Data(),partLabel[1].Data()),"p");
      legendPtResM->AddEntry(h_TMeffi_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]],Form("%s (%s)",labels[3].Data(),partLabel[3].Data()),"p");
      legendPtResM->AddEntry(h_TMeffi_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]],Form("%s (%s)",labels[4].Data(),partLabel[3].Data()),"p");
      legendPtResM->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single particles"),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print(Form("%s/Paper_ClusterTMEfficiencySE_MCE.%s", outputDir.Data(), suffix.Data()));
      
    
    histoDummyEffiMCE->Draw();
    DrawGammaLines(0.15, 70, 1,1,2,kGray+1, 2);
        h_TMeffiCls_recSE_MCE[0][1][maxEtaBinCaloDis[region[0]]]->Draw("same,p");
        h_TMeffiCls_recSE_MCE[1][1][maxEtaBinCaloDis[region[1]]]->Draw("same,p");
        h_TMeffiCls_recSE_MCE[2][1][maxEtaBinCaloDis[region[2]]]->Draw("same,p");
        h_TMeffiCls_recSE_MCE[3][3][maxEtaBinCaloDis[region[3]]]->Draw("same,p");
        h_TMeffiCls_recSE_MCE[4][3][maxEtaBinCaloDis[region[4]]]->Draw("same,p");
      legendPtResM->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single particles"),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{clus}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print(Form("%s/Paper_ClusterTMEfficiencyClusterSE_MCE.%s", outputDir.Data(), suffix.Data()));

    //************************************************************************************************
    //************************ TM effi plot ECals all eta regions ************************************
    //************************************************************************************************
    Double_t arrayBoundariesTMEffEtaX1_4[4];
    Double_t arrayBoundariesTMEffEtaY1_4[3];
    Double_t relMarginginsTMEffEtaX[3];
    Double_t relMarginginsTMEffEtaY[3];
    textSizeLabelsPixel             = 750*0.05;
    
    ReturnCorrectValuesForCanvasScaling(2250,750, 3, 1, 0.04, 0.005, 0.01,0.095,arrayBoundariesTMEffEtaX1_4,arrayBoundariesTMEffEtaY1_4,relMarginginsTMEffEtaX,relMarginginsTMEffEtaY);

    TCanvas* canvasTMEffEtaPaperPlot     = new TCanvas("canvasTMEffEtaPaperPlot","",0,0,2250,750);  // gives the page size
    DrawGammaCanvasSettings( canvasTMEffEtaPaperPlot,  0.0, 0.0, 0.0,0.0);

    TPad* padTMEffEtaEEMC               = new TPad("padTMEffEtaEEMC", "", arrayBoundariesTMEffEtaX1_4[0], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaEEMC, relMarginginsTMEffEtaX[0], relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaY[0], relMarginginsTMEffEtaY[2]);
    padTMEffEtaEEMC->Draw();

    TPad* padTMEffEtaBEMC                = new TPad("padTMEffEtaBEMC", "", arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaBEMC, relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaY[0], relMarginginsTMEffEtaY[2]);
    padTMEffEtaBEMC->Draw();

    TPad* padTMEffEtaFEMC                = new TPad("padTMEffEtaFEMC", "", arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[3], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaFEMC, relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaX[2], relMarginginsTMEffEtaY[0], relMarginginsTMEffEtaY[2]);
    padTMEffEtaFEMC->Draw();

    double textsizeLabelsTMEffEta    = 0;
    double textsizeFacTMEffEta       = 0;
    if (padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) < padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1())){
        textsizeLabelsTMEffEta         = (Double_t)textSizeLabelsPixel/padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) ;
        textsizeFacTMEffEta            = (Double_t)1./padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) ;
    } else {
        textsizeLabelsTMEffEta         = (Double_t)textSizeLabelsPixel/padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1());
        textsizeFacTMEffEta            = (Double_t)1./padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1());
    }
    histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 70,1000,0.0, 1.55);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM} #times A", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.77);
    histoDummyEffiMCE->GetXaxis()->SetNoExponent();
    histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
    histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,70);
    histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);
    TLegend* legendEffiE[10] = {nullptr};
    for (Int_t iSet = 0; iSet < 3; iSet++){
      if(iSet==1){
        padTMEffEtaBEMC->cd();
        padTMEffEtaBEMC->SetLogx();
      } else if(iSet==2){
        padTMEffEtaFEMC->cd();
        padTMEffEtaFEMC->SetLogx();
      } else {
        padTMEffEtaEEMC->cd();
        padTMEffEtaEEMC->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==0)drawLatexAdd( labels[iSet].Data(),0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      else drawLatexAdd( labels[iSet].Data(),0.05,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iSet==0)drawLatexAdd( "a)",0.14,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);
      else if(iSet==1)drawLatexAdd( "b)",0.05,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);
      else drawLatexAdd("c)",0.05,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);

      if(iSet==1){
        double startleg = 0.26;
        drawLatexAdd(labelEnergy,0.95,startleg,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single e^{#pm}"),0.95, startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.95, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,70, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 1;
      Int_t nActiveEta            = maxNEtaBinsFull[region[iSet]]+1;
      legendEffiE[iSet]      = GetAndSetLegend2(0.4, 0.94-(nActiveEta/2*0.85*textSizeLabelsRel), 0.95, 0.94,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      for (Int_t iEta=minEtaBinCaloDis[region[iSet]]; iEta<maxEtaBinCaloDis[region[iSet]]+1;iEta++){
        int iEtaPlot = iEta==maxEtaBinCaloDis[region[iSet]] ? nEta : iEta;
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
        h_TMeffi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");

        if (iEta == maxEtaBinCaloDis[region[iSet]] )
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f", nomEtaRange[iSet][0], nomEtaRange[iSet][1]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasTMEffEtaPaperPlot->SaveAs(Form("%s/TMEffEta_PaperPlot.%s",outputDir.Data(), suffix.Data()));

    //************************************************************************************************
    //************************ TM effi MC plot ECals selected eta regions ************************************
    //************************************************************************************************
    for (Int_t iSet = 0; iSet < 3; iSet++){
      if(iSet==1 ){
        padTMEffEtaBEMC->cd();
        padTMEffEtaBEMC->SetLogx();
      } else if(iSet==2 ){
        padTMEffEtaFEMC->cd();
        padTMEffEtaFEMC->SetLogx();
      } else {
        padTMEffEtaEEMC->cd();
        padTMEffEtaEEMC->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==0)drawLatexAdd( "a)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else if(iSet==1)drawLatexAdd( "b)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else drawLatexAdd("c)",0.94,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);

      if(iSet==2){
        double startleg = 0.92;
        drawLatexAdd(labelEnergy,0.94,startleg,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single e^{#pm}"),0.94, startleg-textSizeLabelsRel*1.05,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.94, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 1;
      Int_t nActiveEta      = 3;
      if (iSet == 0) 
        nActiveEta      = 2;
      double offsetLegend = 0.;
      if (iSet == 0)
        offsetLegend = 0.105;
      
      legendEffiE[iSet]      = GetAndSetLegend2(0.03+offsetLegend, 0.95-((nActiveEta+1)*textSizeLabelsRel), 0.3+offsetLegend, 0.95,textSizeLabelsRel, 0, "", 42, 0.25);
      legendEffiE[iSet]->AddEntry((TObject*)0,labels[iSet].Data()," ");
      for (Int_t iEta=0; iEta<3;iEta++){
        if (!hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta] || (iSet == 0 && iEta == 0)) continue;
        DrawGammaSetMarker(hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta], markerCaloPlotSP[caloID[iSet]][iEta], 1.5*sizeCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta]);
        hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");
        legendEffiE[iSet]->AddEntry(hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEtaSpecCalo[caloID[iSet]][iEta],partEtaSpecCalo[caloID[iSet]][iEta+1]),"p");
      }
      legendEffiE[iSet]->Draw();
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasTMEffEtaPaperPlot->SaveAs(Form("%s/TMEffEta_PaperPlotAlt.%s",outputDir.Data(), suffix.Data()));
    
    //************************************************************************************************
    //************************ TM effi MC plot ECals + HCals selected eta regions ************************************
    //************************************************************************************************
    for (Int_t iSet = 0; iSet < 5; iSet++){
      if(iSet==1 || iSet == 3){
        padTMEffEtaBEMC->cd();
        padTMEffEtaBEMC->SetLogx();
      } else if(iSet==2 || iSet == 4){
        padTMEffEtaFEMC->cd();
        padTMEffEtaFEMC->SetLogx();
      } else {
        padTMEffEtaEEMC->cd();
        padTMEffEtaEEMC->SetLogx();
      }

      if (iSet < 3)histoDummyEffiMCE->DrawCopy();
      if(iSet==0)drawLatexAdd( "a)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else if(iSet==1)drawLatexAdd( "b)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else if(iSet==2)drawLatexAdd("c)",0.94,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);

      if(iSet==0){
        double startleg = 0.92;
        drawLatexAdd(labelEnergy,0.95,startleg,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single e^{#pm} (ECal)/ #pi^{#pm} (HCal)"),0.95, startleg-textSizeLabelsRel*1.05,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.95, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 1;
      if (iSet > 2)
        icurrPID = 3;
      Int_t nActiveEta      = 3;
      if (iSet == 0) 
        nActiveEta      = 2;
      else if (iSet == 3) 
        nActiveEta      = 1;
      double offsetLegend = 0.;
      if (iSet == 0)
        offsetLegend = 0.105;
      else if (iSet == 3 || iSet == 4 )
        offsetLegend = 0.6;
      legendEffiE[iSet]      = GetAndSetLegend2(0.03+offsetLegend, 0.95-((nActiveEta+1)*textSizeLabelsRel), 0.3+offsetLegend, 0.95,textSizeLabelsRel, 0, "", 42, 0.25);
      legendEffiE[iSet]->AddEntry((TObject*)0,labels[iSet].Data()," ");
      for (Int_t iEta=0; iEta<3;iEta++){
        if (!hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta] || (iSet == 0 && iEta == 0)) continue;
        DrawGammaSetMarker(hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta], markerCaloPlotSP[caloID[iSet]][iEta], 1.5*sizeCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta]);
        hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");
        legendEffiE[iSet]->AddEntry(hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEtaSpecCalo[caloID[iSet]][iEta],partEtaSpecCalo[caloID[iSet]][iEta+1]),"p");
      }
      legendEffiE[iSet]->Draw();
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasTMEffEtaPaperPlot->SaveAs(Form("%s/TMEffEta_E+HCal_PaperPlotAlt.%s",outputDir.Data(), suffix.Data()));

    //************************************************************************************************
    //************************ effi MC plot ECals all eta regions ************************************
    //************************************************************************************************
    ReturnCorrectValuesForCanvasScaling(2250,750, 3, 1,0.035, 0.005, 0.01,0.095,arrayBoundariesTMEffEtaX1_4,arrayBoundariesTMEffEtaY1_4,relMarginginsTMEffEtaX,relMarginginsTMEffEtaY);

    TCanvas* canvasTMEffEtaPaperPlot2     = new TCanvas("canvasTMEffEtaPaperPlot","",0,0,2250,750);  // gives the page size
    DrawGammaCanvasSettings( canvasTMEffEtaPaperPlot2,  0., 0., 0.,0.);

    TPad* padTMEffEtaEEMC2               = new TPad("padTMEffEtaEEMC2", "", arrayBoundariesTMEffEtaX1_4[0], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaEEMC2, relMarginginsTMEffEtaX[0], relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaY[0], relMarginginsTMEffEtaY[2]);
    padTMEffEtaEEMC2->Draw();

    TPad* padTMEffEtaBEMC2                = new TPad("padTMEffEtaBEMC2", "", arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaBEMC2, relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaY[0], relMarginginsTMEffEtaY[2]);
    padTMEffEtaBEMC2->Draw();

    TPad* padTMEffEtaFEMC2                = new TPad("padTMEffEtaFEMC2", "", arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[3], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaFEMC2, relMarginginsTMEffEtaX[1], relMarginginsTMEffEtaX[2], relMarginginsTMEffEtaY[0], relMarginginsTMEffEtaY[2]);
    padTMEffEtaFEMC2->Draw();
    histoDummyEffiMCE->GetYaxis()->SetTitleOffset(0.72);
    histoDummyEffiMCE->GetYaxis()->SetTitle("#varepsilon #times A");

    for (Int_t iSet = 0; iSet < 3; iSet++){
      if(iSet==1){
        padTMEffEtaBEMC2->cd();
        padTMEffEtaBEMC2->SetLogx();
      } else if(iSet==2){
        padTMEffEtaFEMC2->cd();
        padTMEffEtaFEMC2->SetLogx();
      } else {
        padTMEffEtaEEMC2->cd();
        padTMEffEtaEEMC2->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==0)drawLatexAdd( labels[iSet].Data(),0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      else drawLatexAdd( labels[iSet].Data(),0.05,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iSet==0)drawLatexAdd( "a)",0.14,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);
      else if(iSet==1)drawLatexAdd( "b)",0.05,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);
      else drawLatexAdd("c)",0.05,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);

      if(iSet==1){
        double startleg = 0.26;
        drawLatexAdd(labelEnergy,0.95,startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single e^{#pm}"),0.95, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
        // drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.95, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 1;
      Int_t nActiveEta            = maxNEtaBinsFull[region[iSet]]+1;
      legendEffiE[iSet]      = GetAndSetLegend2(0.4, 0.94-(nActiveEta/2*0.85*textSizeLabelsRel), 0.95, 0.94,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      for (Int_t iEta=minEtaBinCaloDis[region[iSet]]; iEta<maxEtaBinCaloDis[region[iSet]]+1;iEta++){
        int iEtaPlot = iEta==maxEtaBinCaloDis[region[iSet]] ? nEta : iEta;
        DrawGammaSetMarker(h_effi_recSE_MCE[iSet][icurrPID][iEta], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
        h_effi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");

        if (iEta == maxEtaBinCaloDis[region[iSet]] )
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f", nomEtaRange[iSet][0], nomEtaRange[iSet][1]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasTMEffEtaPaperPlot2->SaveAs(Form("%s/ClusRecEffEta_ECals_PaperPlot.%s",outputDir.Data(), suffix.Data()));

    //************************************************************************************************
    //************************ effi MC plot ECals selected eta regions ************************************
    //************************************************************************************************
    for (Int_t iSet = 0; iSet < 3; iSet++){
      if(iSet==1){
        padTMEffEtaBEMC2->cd();
        padTMEffEtaBEMC2->SetLogx();
      } else if(iSet==2){
        padTMEffEtaFEMC2->cd();
        padTMEffEtaFEMC2->SetLogx();
      } else {
        padTMEffEtaEEMC2->cd();
        padTMEffEtaEEMC2->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==0)drawLatexAdd( "a)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else if(iSet==1)drawLatexAdd( "b)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else drawLatexAdd("c)",0.94,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);

      if(iSet==2){
        double startleg = 0.92;
        drawLatexAdd(labelEnergy,0.94,startleg,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single e^{#pm}"),0.94, startleg-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID        = 1;
      Int_t nActiveEta      = 3;
      double offsetLegend = 0.;
      if (iSet == 0)
        offsetLegend = 0.1;
      
      legendEffiE[iSet]      = GetAndSetLegend2(0.03+offsetLegend, 0.95-((nActiveEta+1)*textSizeLabelsRel), 0.3+offsetLegend, 0.95,textSizeLabelsRel, 0, "", 42, 0.25);
      legendEffiE[iSet]->AddEntry((TObject*)0,labels[iSet].Data()," ");
      for (Int_t iEta=0; iEta<3;iEta++){
        if (!hSP_effi_recSE_MCE[iSet][icurrPID][iEta]) continue;
        DrawGammaSetMarker(hSP_effi_recSE_MCE[iSet][icurrPID][iEta], markerCaloPlotSP[caloID[iSet]][iEta], 1.5*sizeCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta]);
        hSP_effi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");
        legendEffiE[iSet]->AddEntry(hSP_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEtaSpecCalo[caloID[iSet]][iEta],partEtaSpecCalo[caloID[iSet]][iEta+1]),"p");
      }
      legendEffiE[iSet]->Draw();
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasTMEffEtaPaperPlot2->SaveAs(Form("%s/ClusRecEffEta_ECals_PaperPlotAlt.%s",outputDir.Data(), suffix.Data()));

    
    //************************************************************************************************
    //************************ effi MC plot ECals selected eta regions ************************************
    //************************************************************************************************
    for (Int_t iSet = 0; iSet < 5; iSet++){
      if(iSet==1 || iSet==3){
        padTMEffEtaBEMC2->cd();
        padTMEffEtaBEMC2->SetLogx();
      } else if(iSet==2 || iSet==4){
        padTMEffEtaFEMC2->cd();
        padTMEffEtaFEMC2->SetLogx();
      } else {
        padTMEffEtaEEMC2->cd();
        padTMEffEtaEEMC2->SetLogx();
      }

      if(iSet<3) histoDummyEffiMCE->DrawCopy();
      if(iSet==0)drawLatexAdd( "a)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else if(iSet==1)drawLatexAdd( "b)",0.96,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);
      else if (iSet==2) drawLatexAdd("c)",0.94,0.14,textSizeLabelsRel,kFALSE,kFALSE,true);

      if(iSet==0){
        double startleg = 0.92;
        drawLatexAdd(labelEnergy,0.95,startleg,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single e^{#pm} (ECal)/ #pi^{#pm} (HCal)"),0.95, startleg-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID        = 1;
      if (iSet > 2)
        icurrPID        = 3;
      Int_t nActiveEta      = 3;
      if (iSet == 3)
        nActiveEta      = 1;
      double offsetLegend = 0.;
      if (iSet == 0)
        offsetLegend = 0.1;
      if (iSet == 4 || iSet == 3 )
        offsetLegend = 0.6;
      
      legendEffiE[iSet]      = GetAndSetLegend2(0.03+offsetLegend, 0.95-((nActiveEta+1)*textSizeLabelsRel), 0.3+offsetLegend, 0.95,textSizeLabelsRel, 0, "", 42, 0.25);
      legendEffiE[iSet]->AddEntry((TObject*)0,labels[iSet].Data()," ");
      for (Int_t iEta=0; iEta<3;iEta++){
        if (!hSP_effi_recSE_MCE[iSet][icurrPID][iEta]) continue;
        DrawGammaSetMarker(hSP_effi_recSE_MCE[iSet][icurrPID][iEta], markerCaloPlotSP[caloID[iSet]][iEta], 1.5*sizeCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta]);
        hSP_effi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");
        legendEffiE[iSet]->AddEntry(hSP_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEtaSpecCalo[caloID[iSet]][iEta],partEtaSpecCalo[caloID[iSet]][iEta+1]),"p");
      }
      legendEffiE[iSet]->Draw();
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasTMEffEtaPaperPlot2->SaveAs(Form("%s/ClusRecEffEta_E+HCals_PaperPlotAlt.%s",outputDir.Data(), suffix.Data()));
    
    //************************************************************************************************
    //************************ effi MC plot HCals all eta regions ************************************
    //************************************************************************************************
    cout << "moving on to HCAL effi plotting" << endl;
    Double_t arrayBoundariesEffEtaHCalX1_4[4];
    Double_t arrayBoundariesEffEtaHCalY1_4[3];
    Double_t relMarginginsEffEtaHCalX[3];
    Double_t relMarginginsEffEtaHCalY[3];
    textSizeLabelsPixel             = 750*0.05;
    ReturnCorrectValuesForCanvasScaling(1500,750, 2, 1,0.05, 0.005, 0.01,0.095,arrayBoundariesEffEtaHCalX1_4,arrayBoundariesEffEtaHCalY1_4,relMarginginsEffEtaHCalX,relMarginginsEffEtaHCalY);

    TCanvas* canvasEffEtaHCalPaperPlot     = new TCanvas("canvasEffEtaHCalPaperPlot","",0,0,1500,750);  // gives the page size
    DrawGammaCanvasSettings( canvasEffEtaHCalPaperPlot,  0., 0., 0.,0.);
    canvasEffEtaHCalPaperPlot->Draw();
    
    TPad* padEffEtaHCalOHCAL               = new TPad("padEffEtaHCalOHCAL", "", arrayBoundariesEffEtaHCalX1_4[0], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalOHCAL, relMarginginsEffEtaHCalX[0], relMarginginsEffEtaHCalX[1], relMarginginsEffEtaHCalY[0], relMarginginsEffEtaHCalY[2]);
    padEffEtaHCalOHCAL->Draw();

    TPad* padEffEtaHCalLFHCAL                = new TPad("padEffEtaHCalLFHCAL", "", arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[2], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalLFHCAL, relMarginginsEffEtaHCalX[1], relMarginginsEffEtaHCalX[2], relMarginginsEffEtaHCalY[0], relMarginginsEffEtaHCalY[2]);
    padEffEtaHCalLFHCAL->Draw();

    double textsizeLabelsEffEtaHCal    = 0;
    double textsizeFacEffEtaHCal       = 0;
    if (padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) < padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1())){
        textsizeLabelsEffEtaHCal         = (Double_t)textSizeLabelsPixel/padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) ;
        textsizeFacEffEtaHCal            = (Double_t)1./padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) ;
    } else {
        textsizeLabelsEffEtaHCal         = (Double_t)textSizeLabelsPixel/padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1());
        textsizeFacEffEtaHCal            = (Double_t)1./padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1());
    }
    histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 70,1000,0.0, 1.55);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon #times A", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.72);
    histoDummyEffiMCE->GetXaxis()->SetNoExponent();
    histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
    histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);
    for (Int_t iSet = 3; iSet < 5; iSet++){
      if(iSet==3){
        padEffEtaHCalOHCAL->cd();
        padEffEtaHCalOHCAL->SetLogx();
      } else if(iSet==4){
        padEffEtaHCalLFHCAL->cd();
        padEffEtaHCalLFHCAL->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==3)drawLatexAdd( labels[iSet].Data(),0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      else drawLatexAdd( labels[iSet].Data(),0.05,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iSet==3)drawLatexAdd( "d)",0.14,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);
      else drawLatexAdd("e)",0.05,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);

      if(iSet==4){
        double startleg = 0.26;
        drawLatexAdd(labelEnergy,0.94,startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single #pi^{#pm}"),0.94, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,99, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 3;
      Int_t nActiveEta            = maxNEtaBinsFull[region[iSet]]+1;
      legendEffiE[iSet]      = GetAndSetLegend2(0.4, 0.94-(nActiveEta/2*0.85*textSizeLabelsRel), 0.95, 0.94,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      for (Int_t iEta=minEtaBinCaloDis[region[iSet]]; iEta<maxEtaBinCaloDis[region[iSet]]+1;iEta++){
        int iEtaPlot = iEta==maxEtaBinCaloDis[region[iSet]] ? nEta : iEta;
        DrawGammaSetMarker(h_effi_recSE_MCE[iSet][icurrPID][iEta], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
        h_effi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");

        if (iEta == maxEtaBinCaloDis[region[iSet]] )
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f", nomEtaRange[iSet][0], nomEtaRange[iSet][1]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasEffEtaHCalPaperPlot->SaveAs(Form("%s/ClusRecEffEta_HCals_PaperPlot.%s",outputDir.Data(), suffix.Data()));

    
    //************************************************************************************************
    //************************ effi MC plot HCals all eta regions ************************************
    //************************************************************************************************
    for (Int_t iSet = 3; iSet < 5; iSet++){
      if(iSet==3){
        padEffEtaHCalOHCAL->cd();
        padEffEtaHCalOHCAL->SetLogx();
      } else if(iSet==4){
        padEffEtaHCalLFHCAL->cd();
        padEffEtaHCalLFHCAL->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==3)drawLatexAdd( "d)",0.93,0.13,textSizeLabelsRel,kFALSE,kFALSE,false);
      else drawLatexAdd("e)",0.92,0.13,textSizeLabelsRel,kFALSE,kFALSE,false);

      if(iSet==3){
        double startleg = 0.92;
        drawLatexAdd(labelEnergy,0.14,startleg,textSizeLabelsRel,kFALSE,kFALSE,false);
        drawLatexAdd(Form("single #pi^{#pm}"),0.14, startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,false);
      }
      DrawGammaLines(0.15,99, 1., 1., 1, kGray+2, 7);
      Int_t nActiveEta       = 3;
      if (caloID[iSet] == 4 || caloID[iSet] == 5 )
        nActiveEta       = 1;
      Int_t icurrPID = 3;
      double offsetLegend = 0.;
      if (iSet == 4)
        offsetLegend = -0.64;
      legendEffiE[iSet]      = GetAndSetLegend2(0.67+offsetLegend, 0.95-((nActiveEta+1)*textSizeLabelsRel), 0.95+offsetLegend, 0.95,textSizeLabelsRel, 1, "", 42, 0.2);
      legendEffiE[iSet]->AddEntry((TObject*)0,labels[iSet].Data()," ");
      for (Int_t iEta=0; iEta<3;iEta++){
        if (!hSP_effi_recSE_MCE[iSet][icurrPID][iEta]) continue;
        DrawGammaSetMarker(hSP_effi_recSE_MCE[iSet][icurrPID][iEta], markerCaloPlotSP[caloID[iSet]][iEta], 1.5*sizeCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta]);
        hSP_effi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");
        legendEffiE[iSet]->AddEntry(hSP_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEtaSpecCalo[caloID[iSet]][iEta],partEtaSpecCalo[caloID[iSet]][iEta+1]),"p");
      }
      legendEffiE[iSet]->Draw();
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasEffEtaHCalPaperPlot->SaveAs(Form("%s/ClusRecEffEta_HCals_PaperPlotAlt.%s",outputDir.Data(), suffix.Data()));

    //************************************************************************************************
    //************************ TM effi plot HCals all eta regions ************************************
    //************************************************************************************************
    cout << "moving on to HCAL TM effi plotting" << endl;
    textSizeLabelsPixel             = 750*0.05;
    ReturnCorrectValuesForCanvasScaling(1500,750, 2, 1,0.0565, 0.005, 0.01,0.095,arrayBoundariesEffEtaHCalX1_4,arrayBoundariesEffEtaHCalY1_4,relMarginginsEffEtaHCalX,relMarginginsEffEtaHCalY);
    
    TCanvas* canvasEffEtaHCalPaperPlot2     = new TCanvas("canvasEffEtaHCalPaperPlot","",0,0,1500,750);  // gives the page size
    DrawGammaCanvasSettings( canvasEffEtaHCalPaperPlot2,  0.05, 0.01, 0.01,0.095);

    cout << "padEffEtaHCalOHCAL2" << endl;
    TPad* padEffEtaHCalOHCAL2               = new TPad("padEffEtaHCalOHCAL2", "", arrayBoundariesEffEtaHCalX1_4[0], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalOHCAL2, relMarginginsEffEtaHCalX[0], relMarginginsEffEtaHCalX[1], relMarginginsEffEtaHCalY[0], relMarginginsEffEtaHCalY[2]);
    padEffEtaHCalOHCAL2->Draw();

    cout << "padEffEtaHCalLFHCAL2" << endl;
    TPad* padEffEtaHCalLFHCAL2                = new TPad("padEffEtaHCalLFHCAL2", "", arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[2], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalLFHCAL2, relMarginginsEffEtaHCalX[1], relMarginginsEffEtaHCalX[2], relMarginginsEffEtaHCalY[0], relMarginginsEffEtaHCalY[2]);
    padEffEtaHCalLFHCAL2->Draw();

    histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 70,1000,0.0, 1.55);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM} #times A", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.75);
    histoDummyEffiMCE->GetXaxis()->SetNoExponent();
    histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
    histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,50);
    histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);
    for (Int_t iSet = 3; iSet < 5; iSet++){
      if(iSet==3){
        padEffEtaHCalOHCAL2->cd();
        padEffEtaHCalOHCAL2->SetLogx();
      } else {
        padEffEtaHCalLFHCAL2->cd();
        padEffEtaHCalLFHCAL2->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==3)drawLatexAdd( labels[iSet].Data(),0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      else drawLatexAdd( labels[iSet].Data(),0.05,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(iSet==3)drawLatexAdd( "d)",0.14,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);
      else drawLatexAdd("e)",0.05,0.86,textSizeLabelsRel,kFALSE,kFALSE,false);

      if(iSet==4){
        double startleg = 0.26;
        drawLatexAdd(labelEnergy,0.95,startleg,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("single #pi^{#pm}"),0.95, startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
        drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.95, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 3;
      Int_t nActiveEta            = maxNEtaBinsFull[region[iSet]]+1;
      legendEffiE[iSet]      = GetAndSetLegend2(0.4, 0.94-(nActiveEta/2*0.85*textSizeLabelsRel), 0.95, 0.94,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      for (Int_t iEta=minEtaBinCaloDis[region[iSet]]; iEta<maxEtaBinCaloDis[region[iSet]]+1;iEta++){
        int iEtaPlot = iEta==maxEtaBinCaloDis[region[iSet]] ? nEta : iEta;
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
        h_TMeffi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");

        if (iEta == maxEtaBinCaloDis[region[iSet]] )
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f", nomEtaRange[iSet][0], nomEtaRange[iSet][1]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasEffEtaHCalPaperPlot2->SaveAs(Form("%s/TMEffEta_PaperPlot_HCal.%s",outputDir.Data(), suffix.Data()));
    
    //************************************************************************************************
    //************************ TM effi plot HCals selected eta regions ************************************
    //************************************************************************************************
    for (Int_t iSet = 3; iSet < 5; iSet++){
      if(iSet==3){
        padEffEtaHCalOHCAL2->cd();
        padEffEtaHCalOHCAL2->SetLogx();
      } else {
        padEffEtaHCalLFHCAL2->cd();
        padEffEtaHCalLFHCAL2->SetLogx();
      }

      histoDummyEffiMCE->DrawCopy();
      if(iSet==3)drawLatexAdd( "d)",0.93,0.13,textSizeLabelsRel,kFALSE,kFALSE,false);
      else drawLatexAdd("e)",0.92,0.13,textSizeLabelsRel,kFALSE,kFALSE,false);

      if(iSet==3){
        double startleg = 0.92;
        drawLatexAdd(labelEnergy,0.15,startleg,textSizeLabelsRel,kFALSE,kFALSE,false);
        drawLatexAdd(Form("single #pi^{#pm}"),0.15, startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,false);
        drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.15, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,false);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 3;
      Int_t nActiveEta       = 3;
      if (caloID[iSet] == 4 || caloID[iSet] == 5 )
        nActiveEta       = 1;
      double offsetLegend = 0.;
      if (iSet == 4)
        offsetLegend = -0.64;
      legendEffiE[iSet]      = GetAndSetLegend2(0.67+offsetLegend, 0.95-((nActiveEta+1)*textSizeLabelsRel), 0.95+offsetLegend, 0.95,textSizeLabelsRel, 1, "", 42, 0.2);
      legendEffiE[iSet]->AddEntry((TObject*)0,labels[iSet].Data()," ");
      for (Int_t iEta=0; iEta<3;iEta++){
        if (!hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta]) continue;
        DrawGammaSetMarker(hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta], markerCaloPlotSP[caloID[iSet]][iEta], 1.5*sizeCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta], colorCaloPlotSP[caloID[iSet]][iEta]);
        hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");

        legendEffiE[iSet]->AddEntry(hSP_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEtaSpecCalo[caloID[iSet]][iEta],partEtaSpecCalo[caloID[iSet]][iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
    canvasEffEtaHCalPaperPlot2->SaveAs(Form("%s/TMEffEta_PaperPlotAlt_HCal.%s",outputDir.Data(), suffix.Data()));
    
  }
}
