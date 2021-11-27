#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void compare_clustereffi2(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
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
  // for (Int_t iCl = 0; iCl < nClus; iCl++){
  //   gSystem->Exec("mkdir -p "+outputDir+"/"+calo+nameClus);
  // }
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
  Size_t markerSizeSet[8]     = {1.8, 2.0, 2.4, 1.5, 1.8, 1.8, 1.5, 1.5 };
  
  const Int_t maxNSets        = 10;
  TString outTrackCuts[3]     = {"", "LI2", "LI3"};
  TString labelTrackCuts[3]   = {"", "#geq 2 tracker hits", "#geq 3 tracker hits"};

  TString outBetaCuts[2]      = {"", "LI3"};
  TString labelBetaCuts[2]    = {"", "#geq 3 tracker hits"};
  
  
  TH1D* h_effi_rec_E[maxNSets][nPID][nEta+1]     = {{{NULL}}};
  TH1D* h_effi_rec_MCE[maxNSets][nPID][nEta+1]   = {{{NULL}}};
  TH1D* h_effi_recSE_E[maxNSets][nPID][nEta+1]   = {{{NULL}}};
  TH1D* h_effi_recSE_MCE[maxNSets][nPID][nEta+1] = {{{NULL}}};
  TH1D* h_TMeffi_recSE_MCE[maxNSets][nPID][nEta+1] = {{{NULL}}};
  TH1D* h_TMeffiCls_recSE_MCE[maxNSets][nPID][nEta+1] = {{{NULL}}};
  TH1D* h_cluster_NTowerMean_E[maxNSets]         = {NULL};
  TH1D* h_cluster_NClMean_E[maxNSets]            = {NULL};
  TH1D* h_cluster_NClMean_MCE[maxNSets]          = {NULL};
  TFile* inputFiles[maxNSets]                           = {NULL};
  TString inputFilesNames[maxNSets];
  TString labels[maxNSets];
  TString calo[maxNSets];
  int region[maxNSets];
  
  bool currPIDavail[nPID] = {false};

  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> calo[nSets] >> region[nSets]>> labels[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< calo[nSets].Data()<< "\t"<< labels[nSets].Data() <<std::endl;
      labels[nSets].ReplaceAll("_"," ");
      colorSet[nSets] = getCaloColor(calo[nSets]);
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;


  for (Int_t iSet = 0; iSet < nSets; iSet++){
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
        }
      }
      h_cluster_NTowerMean_E[iSet] = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NTowerMean_%s_E", calo[iSet].Data(), "MA"));
      h_cluster_NClMean_E[iSet]    = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NClMean_%s_E", calo[iSet].Data(), "MA"));
      h_cluster_NClMean_MCE[iSet]  = (TH1D*)inputFiles[iSet]->Get(Form("ClusterEfficiencyProjections_%s/h_CS_NClMean_%s_MCE", calo[iSet].Data(), "MA"));
  }

  TCanvas* cReso = new TCanvas("cReso","",0,0,900,800);
  DrawGammaCanvasSettings( cReso, 0.085, 0.025, 0.02, 0.105);
  cReso->SetLogx();

  TH2F* histoDummyEffiE   = new TH2F("histoDummyEffiE","histoDummyEffiE",1000,0.15, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiE, "#it{E} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.87);
  histoDummyEffiE->GetXaxis()->SetNoExponent();
  histoDummyEffiE->GetXaxis()->SetRangeUser(0.15,50);
  histoDummyEffiE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiE->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.57);
  histoDummyEffiMCE->GetXaxis()->SetNoExponent();
  histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,50);
  histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);
  TLegend* legendPtResM  = GetAndSetLegend2(0.14, 0.72-(nSets*textSizeLabelsRel), 0.34, 0.72,textSizeLabelsPixel, 1, "", 43, 0.15);

  for (Int_t pid = 1; pid < nPID; pid++){
    if(!currPIDavail[pid]) continue;
    // histoDummyEffiE->Draw();
      legendPtResM->Clear();
    if(pid==1 && addName.Contains("ECal")) legendPtResM  = GetAndSetLegend2(0.8, 0.935-(nSets*textSizeLabelsRel), 0.94, 0.935,textSizeLabelsPixel, 1, "", 43, 0.25);
    if(pid==1 && addName.Contains("HCal")) legendPtResM  = GetAndSetLegend2(0.75, 0.935-(nSets*textSizeLabelsRel), 0.94, 0.935,textSizeLabelsPixel, 1, "", 43, 0.25);
    
    histoDummyEffiMCE->Draw();
    DrawGammaLines(0.15, 50, 1,1,2,kGray+1, 2);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
        h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]]->Draw("same,p");
        legendPtResM->AddEntry(h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]],labels[iSet].Data(),"p");
      }
      legendPtResM->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single %s", partLabel[pid].Data()),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print(Form("%s/%s_ClusterTMEfficiencySE_MCE.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    
    
    histoDummyEffiMCE->Draw();
    DrawGammaLines(0.15, 50, 1,1,2,kGray+1, 2);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_TMeffiCls_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
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
}
