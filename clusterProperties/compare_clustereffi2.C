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
  Size_t markerSizeSet[8]     = {1.8, 2.0, 2.4, 1.5, 2.4, 1.8, 1.5, 1.5 };
  
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
    // getCaloColor(labels[iSet],false) = getCaloColor(labels[iSet],false);
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
  TH2F* histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 100,1000,0.0, 1.55);
  SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.57);
  histoDummyEffiMCE->GetXaxis()->SetNoExponent();
  histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
  histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,50);
  histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);

  TH2F* histoDummyNTowerMean   = new TH2F("histoDummyNTowerMean","histoDummyNTowerMean",1000,0.21, 99,1000,1, 999);
  SetStyleHistoTH2ForGraphs(histoDummyNTowerMean, "#it{E} (GeV)","#LT#it{N}_{tower}#GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.8);
  histoDummyNTowerMean->GetXaxis()->SetNoExponent();
  histoDummyNTowerMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyNTowerMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyNClPerParMean   = new TH2F("histoDummyNClPerParMean","histoDummyNClPerParMean",1000,0.21, 99,1000,0.0, 4.9);
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
    // histoDummyEffiE->Draw();
      legendPtResM->Clear();
    // if(pid==1 && addName.Contains("ECal")) legendPtResM  = GetAndSetLegend2(0.8, 0.935-(nSets*textSizeLabelsRel), 0.94, 0.935,textSizeLabelsPixel, 1, "", 43, 0.25);
    // if(pid==1 && addName.Contains("HCal")) legendPtResM  = GetAndSetLegend2(0.75, 0.935-(nSets*textSizeLabelsRel), 0.94, 0.935,textSizeLabelsPixel, 1, "", 43, 0.25);
    
    histoDummyNClPerParMean->Draw();
    // DrawGammaLines(0.15, 50, 1,1,2,kGray+1, 2);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_cluster_NClMean_E[iSet], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_cluster_NClMean_E[iSet]->Draw("same,hist,p");
        legendPtResM->AddEntry(h_cluster_NClMean_E[iSet],labels[iSet].Data(),"p");
        legendPtResMLeft->AddEntry(h_cluster_NClMean_E[iSet],labels[iSet].Data(),"p");
      }
      legendPtResM->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      // drawLatexAdd(Form("single %s", partLabel[pid].Data()),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("single particles"),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cReso->Print(Form("%s/%s_NClusterPerParticle_Mean_MCE.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));

    
    
  cReso->SetLogy(true);
    
    histoDummyNTowerMean->Draw();
    // DrawGammaLines(0.15, 50, 1,1,2,kGray+1, 2);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_cluster_NTowerMean_E[iSet], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_cluster_NTowerMean_E[iSet]->Draw("same,hist,p");
        // legendPtResM->AddEntry(h_cluster_NTowerMean_E[iSet],labels[iSet].Data(),"p");
      }
      legendPtResMLeft->Draw();
      drawLatexAdd(labelEnergy,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
      drawLatexAdd(Form("single particles"),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("single %s", partLabel[pid].Data()),0.14, 0.91-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print(Form("%s/%s_NTowerInCluster_Mean_E.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
  cReso->SetLogy(false);
    
    
    histoDummyEffiMCE->Draw();
    DrawGammaLines(0.15, 50, 1,1,2,kGray+1, 2);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
        h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]]->Draw("same,p");
        // legendPtResM->AddEntry(h_TMeffi_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]],labels[iSet].Data(),"p");
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
  bool paperPlot = true;
  if(nSets==5 && paperPlot){





  Double_t arrayBoundariesX1_4[2];
  Double_t arrayBoundariesY1_4[3];
  Double_t relativeMarginsX[3];
  Double_t relativeMarginsY[3];
  textSizeLabelsPixel             = 50;
  ReturnCorrectValuesForCanvasScaling(1350,1450, 1, 2,0.07, 0.006, 0.005,0.075,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

  TCanvas* canvasGeoPaper     = new TCanvas("canvasGeoPaper","",0,0,1350,1450);  // gives the page size
  DrawGammaCanvasSettings( canvasGeoPaper,  0.13, 0.02, 0.03, 0.06);

  TPad* padGeoECal               = new TPad("padGeoECal", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
  DrawGammaPadSettings( padGeoECal, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
  padGeoECal->Draw();

  TPad* padGeoHCal                = new TPad("padGeoHCal", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
  DrawGammaPadSettings( padGeoHCal, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
  padGeoHCal->Draw();

  padGeoECal->cd();
  padGeoECal->SetLogx();

  Double_t margin                 = relativeMarginsX[0]*2.7*1350;
  Double_t textsizeLabelsWidth    = 0;
  Double_t textsizeFacWidth       = 0;
  if (padGeoECal->XtoPixel(padGeoECal->GetX2()) < padGeoECal->YtoPixel(padGeoECal->GetY1())){
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padGeoECal->XtoPixel(padGeoECal->GetX2()) ;
      textsizeFacWidth            = (Double_t)1./padGeoECal->XtoPixel(padGeoECal->GetX2()) ;
  } else {
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padGeoECal->YtoPixel(padGeoECal->GetY1());
      textsizeFacWidth            = (Double_t)1./padGeoECal->YtoPixel(padGeoECal->GetY1());
  }

  TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM",1000,0.21, 99,1000,0.01, 3.9);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{E} (GeV)","#LT#it{N}_{cl}/particle #GT", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.35,512,505);//#it{p}_{T} (GeV/#it{c})
  histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
  // histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
  histo2DAllPi0FWHM->SetMaximum(20);
  histo2DAllPi0FWHM->DrawCopy();

  for(Int_t iSet=0; iSet<nSets;iSet++){
    DrawGammaSetMarker(h_cluster_NClMean_E[iSet], getCaloMarker(labels[iSet]), 1.5*markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
    h_cluster_NClMean_E[iSet]->Draw("same,hist,p");
    // legendPtResM->AddEntry(h_cluster_NClMean_E[iSet],labels[iSet].Data(),"p");
    // legendPtResMLeft2->AddEntry(h_cluster_NClMean_E[iSet],labels[iSet].Data(),"p");
  }
  // legendPhiEta    = GetAndSetLegend2(0.10, 0.08, 0.38, 0.08+(3*0.85*textsizeLabelsWidth),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  // legendPhiEta->AddEntry(h_PhiEta_ECal, "all towers","plf");
  // legendPhiEta->AddEntry(h_PhiEta_iEta_ECal, "towers iEta = 1, 30, 70","fl");
  // legendPhiEta->AddEntry(h_PhiEta_iPhi_ECal, "towers iPhi = 1, 50, 100","fl");
  // legendPhiEta->Draw();
  
  drawLatexAdd("#it{#bf{ECCE}} simulation",0.10,0.86,textsizeLabelsWidth,kFALSE,kFALSE,false);
  drawLatexAdd(Form("single particles"),0.10,0.78,textsizeLabelsWidth,kFALSE,kFALSE,false);
  drawLatexAdd(Form("a)"),0.97,0.08,textsizeLabelsWidth,kFALSE,kFALSE,true);
  padGeoHCal->SetLogx();
  padGeoHCal->SetLogy();

  padGeoHCal->cd();
  Double_t textsizeLabelsMass         = 0;
  Double_t textsizeFacMass            = 0;
  if (padGeoHCal->XtoPixel(padGeoHCal->GetX2()) <padGeoHCal->YtoPixel(padGeoHCal->GetY1()) ){
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padGeoHCal->XtoPixel(padGeoHCal->GetX2()) ;
      textsizeFacMass                 = (Double_t)1./padGeoHCal->XtoPixel(padGeoHCal->GetX2()) ;
  } else {
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padGeoHCal->YtoPixel(padGeoHCal->GetY1());
      textsizeFacMass                 = (Double_t)1./padGeoHCal->YtoPixel(padGeoHCal->GetY1());
  }

  TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",1000,0.21, 99,1000,1, 999);//, 100.1, 160.9);//125.1, 155.9);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{E} (GeV)","#LT#it{N}_{tower}#GT", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.35,512,505);
  histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0Mass->GetXaxis()->SetNoExponent();
  histo2DAllPi0Mass->SetMaximum(20);
  histo2DAllPi0Mass->DrawCopy();
  TLegend* legendPtResMLeft2  = GetAndSetLegend2(0.10, 0.95-(nSets*textsizeLabelsMass), 0.28, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);

  for(Int_t iSet=0; iSet<nSets;iSet++){
    DrawGammaSetMarker(h_cluster_NTowerMean_E[iSet], getCaloMarker(labels[iSet]), 1.5*markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
    h_cluster_NTowerMean_E[iSet]->Draw("same,hist,p");
    legendPtResMLeft2->AddEntry(h_cluster_NTowerMean_E[iSet],labels[iSet].Data(),"p");
  }
   legendPtResMLeft2->Draw();
 // drawLatexAdd(Form("ECals"),0.10,0.88,textsizeLabelsMass,kFALSE,kFALSE,false);
  drawLatexAdd(Form("b)"),0.97,0.19,textsizeLabelsMass,kFALSE,kFALSE,true);

  canvasGeoPaper->Update();
  canvasGeoPaper->Print(Form("%s/NTowAndMeanClus_Paper.%s",outputDir.Data(),suffix.Data()));


  padGeoECal->cd();
  padGeoECal->SetLogx();

  margin                 = relativeMarginsX[0]*2.7*1350;
  textsizeLabelsWidth    = 0;
  textsizeFacWidth       = 0;
  if (padGeoECal->XtoPixel(padGeoECal->GetX2()) < padGeoECal->YtoPixel(padGeoECal->GetY1())){
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padGeoECal->XtoPixel(padGeoECal->GetX2()) ;
      textsizeFacWidth            = (Double_t)1./padGeoECal->XtoPixel(padGeoECal->GetX2()) ;
  } else {
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padGeoECal->YtoPixel(padGeoECal->GetY1());
      textsizeFacWidth            = (Double_t)1./padGeoECal->YtoPixel(padGeoECal->GetY1());
  }

  TH2F * histo2DAllEffiTM2    = new TH2F("histo2DAllEffiTM2","histo2DAllEffiTM2",1000,0.26, 55,1000,0.001, 1.35);
  SetStyleHistoTH2ForGraphs(histo2DAllEffiTM2, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.45,512,505);//#it{p}_{T} (GeV/#it{c})
  histo2DAllEffiTM2->GetYaxis()->SetMoreLogLabels(kTRUE);
  // histo2DAllEffiTM2->GetYaxis()->SetNdivisions(505);
  histo2DAllEffiTM2->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllEffiTM2->GetXaxis()->SetTickLength(0.05);
  histo2DAllEffiTM2->GetYaxis()->SetTickLength(0.026);
  histo2DAllEffiTM2->SetMaximum(20);
  histo2DAllEffiTM2->DrawCopy();
  DrawGammaLines(0.26, 55, 1,1,2,kGray+1, 2);
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

  // drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  
  drawLatexAdd("#it{#bf{ECCE}} simulation",0.10,0.86,textsizeLabelsWidth,kFALSE,kFALSE,false);
  drawLatexAdd(Form("single particles"),0.10,0.78,textsizeLabelsWidth,kFALSE,kFALSE,false);
  drawLatexAdd(Form("a)"),0.97,0.08,textsizeLabelsWidth,kFALSE,kFALSE,true);
  histo2DAllEffiTM2->Draw("same,axis");
  padGeoHCal->SetLogx();
  padGeoHCal->SetLogy(false);

  padGeoHCal->cd();
  textsizeLabelsMass         = 0;
  textsizeFacMass            = 0;
  if (padGeoHCal->XtoPixel(padGeoHCal->GetX2()) <padGeoHCal->YtoPixel(padGeoHCal->GetY1()) ){
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padGeoHCal->XtoPixel(padGeoHCal->GetX2()) ;
      textsizeFacMass                 = (Double_t)1./padGeoHCal->XtoPixel(padGeoHCal->GetX2()) ;
  } else {
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padGeoHCal->YtoPixel(padGeoHCal->GetY1());
      textsizeFacMass                 = (Double_t)1./padGeoHCal->YtoPixel(padGeoHCal->GetY1());
  }

  TH2F * histo2DAllEffiTM3            = new TH2F("histo2DAllEffiTM3","histo2DAllEffiTM3",1000,0.26, 55,1000,0.0, 1.35);//, 100.1, 160.9);//125.1, 155.9);
  SetStyleHistoTH2ForGraphs(histo2DAllEffiTM3, "#it{E}^{MC} (GeV)","#kappa_{TM}", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.45,512,505);
  histo2DAllEffiTM3->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllEffiTM3->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllEffiTM3->GetXaxis()->SetTickLength(0.05);
  histo2DAllEffiTM3->GetXaxis()->SetNoExponent();
  histo2DAllEffiTM3->SetMaximum(20);
  histo2DAllEffiTM3->DrawCopy();
  legendPtResMLeft2  = GetAndSetLegend2(0.10, 0.94-(2*textsizeLabelsMass), 0.78, 0.94,textSizeLabelsPixel, 3, "", 43, 0.15);
  for (Int_t pid = 1; pid < nPID; pid++){
    if(!currPIDavail[pid]) continue;
    for(Int_t iSet=0; iSet<nSets;iSet++){
        DrawGammaSetMarker(h_TMeffiCls_recSE_MCE[iSet][pid][maxEtaBinCaloDis[region[iSet]]], getCaloMarker(labels[iSet]), 1.5*markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
    }
  }
    DrawGammaLines(0.26, 55, 1,1,2,kGray+1, 2);
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
      // drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{clus}"),0.14, 0.91-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  
  // for(Int_t iSet=0; iSet<nSets;iSet++){
  //   DrawGammaSetMarker(h_cluster_NTowerMean_E[iSet], getCaloMarker(labels[iSet]), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
  //   h_cluster_NTowerMean_E[iSet]->Draw("same,hist,p");
  //   legendPtResMLeft2->AddEntry(h_cluster_NTowerMean_E[iSet],labels[iSet].Data(),"p");
  // }
   legendPtResMLeft2->Draw();
 // drawLatexAdd(Form("ECals"),0.10,0.88,textsizeLabelsMass,kFALSE,kFALSE,false);
  drawLatexAdd(Form("b)"),0.97,0.19,textsizeLabelsMass,kFALSE,kFALSE,true);
  histo2DAllEffiTM3->Draw("same,axis");

  canvasGeoPaper->Update();
  canvasGeoPaper->Print(Form("%s/ClusterTMEfficiency_Paper.%s",outputDir.Data(),suffix.Data()));








  histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.55);
    histoDummyEffiMCE->Draw();
    legendPtResM  = GetAndSetLegend2(0.68, 0.935-(nSets*textSizeLabelsRel), 0.86, 0.935,textSizeLabelsPixel, 1, "", 43, 0.25);
    legendPtResM->SetMargin(0.3);
    DrawGammaLines(0.15, 50, 1,1,2,kGray+1, 2);
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
    DrawGammaLines(0.15, 50, 1,1,2,kGray+1, 2);
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



    Double_t arrayBoundariesTMEffEtaX1_4[4];
    Double_t arrayBoundariesTMEffEtaY1_4[3];
    Double_t relativeMarginsTMEffEtaX[3];
    Double_t relativeMarginsTMEffEtaY[3];
    textSizeLabelsPixel             = 750*0.05;
    ReturnCorrectValuesForCanvasScaling(2250,750, 3, 1,0.035, 0.01, 0.01,0.095,arrayBoundariesTMEffEtaX1_4,arrayBoundariesTMEffEtaY1_4,relativeMarginsTMEffEtaX,relativeMarginsTMEffEtaY);

    TCanvas* canvasTMEffEtaPaperPlot     = new TCanvas("canvasTMEffEtaPaperPlot","",0,0,2250,750);  // gives the page size
    DrawGammaCanvasSettings( canvasTMEffEtaPaperPlot,  0.05, 0.01, 0.01,0.095);
    // canvasTMEffEtaPaperPlot->SetLogy(1);

    TPad* padTMEffEtaEEMC               = new TPad("padTMEffEtaEEMC", "", arrayBoundariesTMEffEtaX1_4[0], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaEEMC, relativeMarginsTMEffEtaX[0], relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaY[0], relativeMarginsTMEffEtaY[2]);
    padTMEffEtaEEMC->Draw();

    TPad* padTMEffEtaBEMC                = new TPad("padTMEffEtaBEMC", "", arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaBEMC, relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaY[0], relativeMarginsTMEffEtaY[2]);
    padTMEffEtaBEMC->Draw();

    TPad* padTMEffEtaFEMC                = new TPad("padTMEffEtaFEMC", "", arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[3], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaFEMC, relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaX[2], relativeMarginsTMEffEtaY[0], relativeMarginsTMEffEtaY[2]);
    padTMEffEtaFEMC->Draw();

    // padTMEffEtaEEMC->cd();
    // padTMEffEtaEEMC->SetLogy();

    margin                 = relativeMarginsTMEffEtaX[0]*2.7*1350;
    double textsizeLabelsTMEffEta    = 0;
    double textsizeFacTMEffEta       = 0;
    if (padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) < padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1())){
        textsizeLabelsTMEffEta         = (Double_t)textSizeLabelsPixel/padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) ;
        textsizeFacTMEffEta            = (Double_t)1./padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) ;
    } else {
        textsizeLabelsTMEffEta         = (Double_t)textSizeLabelsPixel/padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1());
        textsizeFacTMEffEta            = (Double_t)1./padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1());
    }
    histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 100,1000,0.0, 1.55);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.57);
    histoDummyEffiMCE->GetXaxis()->SetNoExponent();
    histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
    histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,50);
    histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);
    // histo2DPi0TMEffEtaDummy->GetYaxis()->SetTitle("norm. counts");
    // histo2DPi0TMEffEtaDummy->GetYaxis()->SetTitleOffset(1.0);
    // histo2DPi0TMEffEtaDummy->GetXaxis()->SetTickLength(0.025);
    // histo2DPi0TMEffEtaDummy->GetXaxis()->SetRangeUser(0,1.99);
    // histo2DPi0TMEffEtaDummy->GetYaxis()->SetRangeUser(0.1,2.9);
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
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      Int_t icurrPID = 1;
      Int_t nActiveEta            = maxNEtaBinsFull[region[iSet]]+1;
      legendEffiE[iSet]      = GetAndSetLegend2(0.4, 0.94-(nActiveEta/2*0.85*textSizeLabelsRel), 0.95, 0.94,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      for (Int_t iEta=minEtaBinCaloDis[region[iSet]]; iEta<maxEtaBinCaloDis[region[iSet]]+1;iEta++){
        int iEtaPlot = iEta==maxEtaBinCaloDis[region[iSet]] ? nEta : iEta;
        DrawGammaSetMarker(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
        h_TMeffi_recSE_MCE[iSet][icurrPID][iEta]->Draw("same");

        if (iEta == maxEtaBinCaloDis[region[iSet]] )
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[minEtaBinCaloDis[region[iSet]]],partEta[maxEtaBinCaloDis[region[iSet]]]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
      canvasTMEffEtaPaperPlot->SaveAs(Form("%s/TMEffEta_PaperPlot.%s",outputDir.Data(), suffix.Data()));


    histoDummyEffiMCE->GetYaxis()->SetTitle("#varepsilon");

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
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[minEtaBinCaloDis[region[iSet]]],partEta[maxEtaBinCaloDis[region[iSet]]]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
      canvasTMEffEtaPaperPlot->SaveAs(Form("%s/ClusRecEffEta_ECals_PaperPlot.%s",outputDir.Data(), suffix.Data()));





    Double_t arrayBoundariesEffEtaHCalX1_4[4];
    Double_t arrayBoundariesEffEtaHCalY1_4[3];
    Double_t relativeMarginsEffEtaHCalX[3];
    Double_t relativeMarginsEffEtaHCalY[3];
    textSizeLabelsPixel             = 750*0.05;
    ReturnCorrectValuesForCanvasScaling(1500,750, 2, 1,0.045, 0.01, 0.01,0.095,arrayBoundariesEffEtaHCalX1_4,arrayBoundariesEffEtaHCalY1_4,relativeMarginsEffEtaHCalX,relativeMarginsEffEtaHCalY);

    TCanvas* canvasEffEtaHCalPaperPlot     = new TCanvas("canvasEffEtaHCalPaperPlot","",0,0,1500,750);  // gives the page size
    DrawGammaCanvasSettings( canvasEffEtaHCalPaperPlot,  0.05, 0.01, 0.01,0.095);
    // canvasEffEtaHCalPaperPlot->SetLogy(1);

    TPad* padEffEtaHCalOHCAL               = new TPad("padEffEtaHCalOHCAL", "", arrayBoundariesEffEtaHCalX1_4[0], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalOHCAL, relativeMarginsEffEtaHCalX[0], relativeMarginsEffEtaHCalX[1], relativeMarginsEffEtaHCalY[0], relativeMarginsEffEtaHCalY[2]);
    padEffEtaHCalOHCAL->Draw();

    TPad* padEffEtaHCalLFHCAL                = new TPad("padEffEtaHCalLFHCAL", "", arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[2], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalLFHCAL, relativeMarginsEffEtaHCalX[1], relativeMarginsEffEtaHCalX[2], relativeMarginsEffEtaHCalY[0], relativeMarginsEffEtaHCalY[2]);
    padEffEtaHCalLFHCAL->Draw();

    // TPad* padEffEtaHCalFEMC                = new TPad("padEffEtaHCalFEMC", "", arrayBoundariesEffEtaHCalX1_4[2], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[3], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    // DrawGammaPadSettings( padEffEtaHCalFEMC, relativeMarginsEffEtaHCalX[1], relativeMarginsEffEtaHCalX[2], relativeMarginsEffEtaHCalY[0], relativeMarginsEffEtaHCalY[2]);
    // padEffEtaHCalFEMC->Draw();

    // padEffEtaHCalOHCAL->cd();
    // padEffEtaHCalOHCAL->SetLogy();

    // Double_t margin                 = relativeMarginsEffEtaHCalX[0]*2.7*1350;
    double textsizeLabelsEffEtaHCal    = 0;
    double textsizeFacEffEtaHCal       = 0;
    if (padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) < padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1())){
        textsizeLabelsEffEtaHCal         = (Double_t)textSizeLabelsPixel/padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) ;
        textsizeFacEffEtaHCal            = (Double_t)1./padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) ;
    } else {
        textsizeLabelsEffEtaHCal         = (Double_t)textSizeLabelsPixel/padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1());
        textsizeFacEffEtaHCal            = (Double_t)1./padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1());
    }


    histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 99,1000,0.0, 1.55);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.57);
    histoDummyEffiMCE->GetXaxis()->SetNoExponent();
    histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
    // histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,50);
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
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[minEtaBinCaloDis[region[iSet]]],partEta[maxEtaBinCaloDis[region[iSet]]]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_effi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
      canvasEffEtaHCalPaperPlot->SaveAs(Form("%s/ClusRecEffEta_HCals_PaperPlot.%s",outputDir.Data(), suffix.Data()));


    histoDummyEffiMCE   = new TH2F("histoDummyEffiMCE","histoDummyEffiMCE",1000,0.15, 100,1000,0.0, 1.55);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCE, "#it{E}^{MC} (GeV)","#varepsilon_{TM}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,1.3*textSizeSinglePad, 0.9,0.57);
    histoDummyEffiMCE->GetXaxis()->SetNoExponent();
    histoDummyEffiMCE->GetYaxis()->SetRangeUser(0.0,1.35);
    histoDummyEffiMCE->GetXaxis()->SetRangeUser(0.15,50);
    histoDummyEffiMCE->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiMCE->GetXaxis()->SetMoreLogLabels(kTRUE);
    // histo2DPi0TMEffEtaDummy->GetYaxis()->SetTitle("norm. counts");
    // histo2DPi0TMEffEtaDummy->GetYaxis()->SetTitleOffset(1.0);
    // histo2DPi0TMEffEtaDummy->GetXaxis()->SetTickLength(0.025);
    // histo2DPi0TMEffEtaDummy->GetXaxis()->SetRangeUser(0,1.99);
    // histo2DPi0TMEffEtaDummy->GetYaxis()->SetRangeUser(0.1,2.9);
    for (Int_t iSet = 3; iSet < 5; iSet++){
      if(iSet==3){
        padEffEtaHCalOHCAL->cd();
        padEffEtaHCalOHCAL->SetLogx();
      } else {
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
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[minEtaBinCaloDis[region[iSet]]],partEta[maxEtaBinCaloDis[region[iSet]]]),"p");
        else 
          legendEffiE[iSet]->AddEntry(h_TMeffi_recSE_MCE[iSet][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
        legendEffiE[iSet]->Draw();
      }
      histoDummyEffiMCE->Draw("same,axis");
    }
      canvasEffEtaHCalPaperPlot->SaveAs(Form("%s/TMEffEta_PaperPlot_HCal.%s",outputDir.Data(), suffix.Data()));

  }
}
