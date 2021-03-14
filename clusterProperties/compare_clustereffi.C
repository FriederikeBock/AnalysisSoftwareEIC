#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void compare_clustereffi(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString calo              = "",
                            TString suffix            = "pdf"
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();


  TString detLabel = "";

//   Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;
  Double_t maxBetaSigma       = 0.05;

  TString outputDir                 = Form("plots/%s/Compare%s",dateForOutput.Data(), addName.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  for (Int_t iCl = 0; iCl < nClus; iCl++){
    gSystem->Exec("mkdir -p "+outputDir+"/"+calo+nameClus[iCl]);
  }
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

  Color_t colorSet[8]         = {kGray+1, kBlack, kRed-6, kRed+2, kGreen+1, kCyan+1, kAzure+2, kBlack };
  Style_t markerStyleSet[8]   = {24, 20, 25,  21, 30, 42, 46, 20};
  Size_t markerSizeSet[8]     = {1.5, 1.4, 1.6, 1.5, 1.8, 1.8, 1.5, 1.5 };
  
  const Int_t maxNSets        = 10;
  TString outTrackCuts[3]     = {"", "LI2", "LI3"};
  TString labelTrackCuts[3]   = {"", "#geq 2 tracker hits", "#geq 3 tracker hits"};

  TString outBetaCuts[2]      = {"", "LI3"};
  TString labelBetaCuts[2]    = {"", "#geq 3 tracker hits"};
  
  
  TH1D* h_effi_rec_E[maxNSets][nPID][nEta+1][nClus]     = {{{{NULL}}}};
  TH1D* h_effi_rec_MCE[maxNSets][nPID][nEta+1][nClus]   = {{{{NULL}}}};
  TH1D* h_effi_recSE_E[maxNSets][nPID][nEta+1][nClus]   = {{{{NULL}}}};
  TH1D* h_effi_recSE_MCE[maxNSets][nPID][nEta+1][nClus] = {{{{NULL}}}};
  TH1D* h_cluster_NTowerMean_E[maxNSets][nClus]         = {{NULL}};
  TH1D* h_cluster_NClMean_E[maxNSets][nClus]            = {{NULL}};
  TH1D* h_cluster_NClMean_MCE[maxNSets][nClus]          = {{NULL}};
  TFile* inputFiles[maxNSets]                           = {NULL};
  TString inputFilesNames[maxNSets];
  TString labels[maxNSets];
  
  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> labels[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< labels[nSets].Data() <<std::endl;
      labels[nSets].ReplaceAll("_"," ");
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;


  for (Int_t iSet = 0; iSet < nSets; iSet++){
    inputFiles[iSet]  = new TFile(inputFilesNames[iSet].Data());
    for (Int_t iCl = 0; iCl < nClus; iCl++){
      for (Int_t iEta = 0; iEta < nEta+1; iEta++){
        for (Int_t pid = 0; pid < nPID; pid++){
            h_effi_rec_E[iSet][pid][iEta][iCl]      = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effi%s_E_%d_%s", calo.Data(), partName[pid].Data(), iEta, nameClus[iCl].Data()));
            h_effi_rec_MCE[iSet][pid][iEta][iCl]    = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effi%s_MCE_%d_%s", calo.Data(), partName[pid].Data(), iEta, nameClus[iCl].Data()));
            h_effi_recSE_E[iSet][pid][iEta][iCl]    = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiSE%s_E_%d_%s", calo.Data(), partName[pid].Data(), iEta, nameClus[iCl].Data()));
            h_effi_recSE_MCE[iSet][pid][iEta][iCl]  = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/effiSE%s_MCE_%d_%s", calo.Data(), partName[pid].Data(), iEta, nameClus[iCl].Data()));
        }
      }
      h_cluster_NTowerMean_E[iSet][iCl] = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NTowerMean_%s_E", calo.Data(), nameClus[iCl].Data());
      h_cluster_NTowerMean_E[iSet][iCl] = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NClMean_%s_E", calo.Data(), nameClus[iCl].Data());
      h_cluster_NTowerMean_E[iSet][iCl] = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NClMean_%s_MCE", calo.Data(), nameClus[iCl].Data());
    }
  }
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    // 1D PLOT
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }

    TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
    DrawGammaCanvasSettings( cReso, 0.085, 0.025, 0.02, 0.105);
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
    
    TLegend* legendPtResM  = GetAndSetLegend2(0.12, 0.94-(nSets/2*0.75*textSizeLabelsRel), 0.4, 0.94,0.75*textSizeLabelsPixel, 1, "", 43, 0.15);
    
  
    for (Int_t iCl = 0; iCl < nClus; iCl++){
      for (Int_t pid = 0; pid < nPID; pid++){
        histoDummyEffiE->Draw();
          legendPtResM->Clear();
          for(Int_t iSet=0; iSet<nSets;iSet++){
            DrawGammaSetMarker(h_effi_rec_E[iSet][pid][iEta][iCl], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
            h_effi_rec_E[iSet][pid][iEta][iCl]->Draw("same,p");
            legendPtResM->AddEntry(h_effi_rec_E[iSet][pid][iEta][iCl],labels[iSet].Data(),"p");
          }
          legendPtResM->Draw();
          drawLatexAdd(collisionSystem,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
          drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%s, %s clusterizer", calo.Data(), nameClus[iCl].Data()),0.14, 0.91-2*nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        cReso->Print(Form("%s/%s%s/%s_ClusterEfficiency_E_%d_%d.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
        
        histoDummyEffiMCE->Draw();
          for(Int_t iSet=0; iSet<nSets;iSet++){
            DrawGammaSetMarker(h_effi_rec_MCE[iSet][pid][iEta][iCl], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
            h_effi_rec_MCE[iSet][pid][iEta][iCl]->Draw("same,p");
          }
          legendPtResM->Draw();
          drawLatexAdd(collisionSystem,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
          drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%s, %s clusterizer", calo.Data(), nameClus[iCl].Data()),0.14, 0.91-2*nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        cReso->Print(Form("%s/%s%s/%s_ClusterEfficiency_E_%d_%d.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10),
        
        histoDummyEffiE->Draw();
          for(Int_t iSet=0; iSet<nSets;iSet++){
            DrawGammaSetMarker(h_effi_recSE_E[iSet][pid][iEta][iCl], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
            h_effi_recSE_E[iSet][pid][iEta][iCl]->Draw("same,p");
          }
          legendPtResM->Draw();
          drawLatexAdd(collisionSystem,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
          drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%s, %s clusterizer, single entry", calo.Data(), nameClus[iCl].Data()),0.14, 0.91-2*nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        cReso->Print(Form("%s/%s%s/%s_ClusterEfficiencySE_E_%d_%d.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
        
        histoDummyEffiMCE->Draw();
          for(Int_t iSet=0; iSet<nSets;iSet++){
            DrawGammaSetMarker(h_effi_recSE_MCE[iSet][pid][iEta][iCl], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
            h_effi_recSE_MCE[iSet][pid][iEta][iCl]->Draw("same,p");
          }
          legendPtResM->Draw();
          drawLatexAdd(collisionSystem,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
          drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%s, %s clusterizer, single entry", calo.Data(), nameClus[iCl].Data()),0.14, 0.91-2*nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        cReso->Print(Form("%s/%s%s/%s_ClusterEfficiencySE_E_%d_%d.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10),
        
      }
    }
  }
  
  TH2F* histoDummyNTowerMean   = new TH2F("histoDummyNTowerMean","histoDummyNTowerMean",1000,0, 200,1000,0.0, 40);
  SetStyleHistoTH2ForGraphs(histoDummyNTowerMean, "#it{E} (GeV)","#LT#it{N}_{tower}#GT", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNTowerMean->GetXaxis()->SetNoExponent();
  histoDummyNTowerMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNTowerMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyNClPerParMean   = new TH2F("histoDummyNClPerParMean","histoDummyNClPerParMean",1000,0, 200,1000,0.0, 15);
  SetStyleHistoTH2ForGraphs(histoDummyNClPerParMean, "#it{E} (GeV)","#LT#it{N}_{cl}/particle #GT", 0.85*textSizeSinglePad,textSizeSinglePad, textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNClPerParMean->GetXaxis()->SetNoExponent();
  histoDummyNClPerParMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNClPerParMean->GetXaxis()->SetMoreLogLabels(kTRUE);

  for (Int_t iCl = 0; iCl < nClus; iCl++){

    TLegend* legendNTowerCl   = GetAndSetLegend2(0.12,  0.94-(nActiceCl/2*0.85*textSizeLabelsRel), 0.4,  0.94, 0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
    if (calo.CompareTo("FHCAL") == 0) histoDummyNTowerMean->GetYaxis()->SetRangeUser(0, 30);
    histoDummyNTowerMean->Draw();
    for (Int_t iSet = 0; iSet < nClus; iSet++){
      DrawGammaSetMarker(h_cluster_NTowerMean_E[iSet][iCl], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      h_cluster_NTowerMean_E[iSet][iCl]->Draw("same,p");      
      legendNTowerCl->AddEntry(h_cluster_NTowerMean_E[iSet][iCl],nameClus[iCl].Data(),"p");
    }
    legendNTowerCl->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s, %s clusterizer, single entry", calo.Data(), nameClus[iCl].Data()),0.14, 0.91-2*nLinesCol*textSizeLabelsRel*1.1,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cReso->Print(Form("%s/%s%s/NTowerInCluster_Mean_E.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), suffix.Data()));

    
    if (calo.CompareTo("FHCAL") == 0) histoDummyNClPerParMean->GetYaxis()->SetRangeUser(0, 10);
    histoDummyNClPerParMean->Draw();
    DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
    for (Int_t iSet = 0; iSet < nClus; iSet++){
      DrawGammaSetMarker(h_cluster_NClMean_E[iSet][iCl], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      h_cluster_NClMean_E[iSet][iCl]->Draw("same,p");      
    }
    legendNTowerCl->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s, %s clusterizer, single entry", calo.Data(), nameClus[iCl].Data()),0.14, 0.91-2*nLinesCol*textSizeLabelsRel*1.1,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);  
    cReso->Print(Form("%s/%s%s/NClusterPerParticle_Mean_E.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), suffix.Data()));

    histoDummyNClPerParMean->SetXTitle("#it{E}_{MC} (GeV)");
    histoDummyNClPerParMean->Draw();
    DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
    for (Int_t iSet = 0; iSet < nClus; iSet++){
      DrawGammaSetMarker(h_cluster_NClMean_MCE[iSet][iCl], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      h_cluster_NClMean_MCE[iSet][iCl]->Draw("same,p");      
    }
    legendNTowerCl->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s, %s clusterizer, single entry", calo.Data(), nameClus[iCl].Data()),0.14, 0.91-2*nLinesCol*textSizeLabelsRel*1.1,0.85*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    
    cReso->Print(Form("%s/%s%s/NClusterPerParticle_Mean_MCE.%s", outputDir.Data(), calo.Data(), nameClus[iCl].Data(), suffix.Data()));
  }
}
