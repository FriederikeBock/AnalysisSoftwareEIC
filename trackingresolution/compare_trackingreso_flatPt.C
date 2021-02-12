#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void compare_trackingreso_flatPt(
                            TString configInputFiles  = "file.txt",
                            TString suffix            = "pdf",
                            Bool_t properFit          = kTRUE
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();


  TString detLabel = "";
  Bool_t enablePlot[8]        = {1, 1, 1, 1, 1, 1, 1, 1};
  Int_t nActiveEta            = 8;
  Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;
  TString addSigma            = "FitSigma";
  TString addMean             = "FitMean";
  TString addOut              = "Fit";
  if (!properFit){
    maxPtSigma          = 0.42;
    maxEtaSigma         = 0.01;
    maxPhiSigma         = 0.1;
    addSigma            = "sigma";
    addMean             = "mean";
    addOut              = "Hist";
  }

  TString outputDir                 = Form("plots/%s/Compare%s",dateForOutput.Data(), addOut.Data());
  gSystem->Exec("mkdir -p "+outputDir);


  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  const Int_t nEta                  = 7;
  Double_t partEta[8]               = { 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};                                    

  Color_t colorEta[8]         = {kRed+1, kOrange+7, kOrange, kSpring+5, kGreen+1, kCyan+1, kAzure+2, kBlack };
  Style_t markerStyleEta[8]   = {24, 25, 27, 28, 30, 42, 46, 20};
  Size_t markerSizeEta[8]     = {1.5, 1.4, 1.9, 1.5, 1.8, 1.8, 1.5, 1.5 };

  Color_t colorSet[8]         = {kGray+1, kBlack, kRed-6, kRed+2, kGreen+1, kCyan+1, kAzure+2, kBlack };
  Style_t markerStyleSet[8]   = {24, 20, 25,  21, 30, 42, 46, 20};
  Size_t markerSizeSet[8]     = {1.5, 1.4, 1.6, 1.5, 1.8, 1.8, 1.5, 1.5 };
  
  const Int_t maxNSets        = 10;
  
  TH1D* h_tracks_mean_pt_reso[maxNSets][nEta+1]       = {{NULL}};
  TH1D* h_tracks_sigma_pt_reso[maxNSets][nEta+1]      = {{NULL}};
  TH1D* h_tracks_mean_pt_resoEta[maxNSets][nEta+1]    = {{NULL}};
  TH1D* h_tracks_sigma_pt_resoEta[maxNSets][nEta+1]   = {{NULL}};
  TH1D* h_tracks_mean_pt_resoPhi[maxNSets][nEta+1]    = {{NULL}};
  TH1D* h_tracks_sigma_pt_resoPhi[maxNSets][nEta+1]   = {{NULL}};
  TFile* inputFiles[maxNSets]                 = {NULL};
  TString inputFilesNames[maxNSets];
  TString labels[maxNSets];
  
  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> labels[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< labels[nSets].Data() <<std::endl;
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;

  for (Int_t iSet = 0; iSet < nSets; iSet++){
    inputFiles[iSet]  = new TFile(inputFilesNames[iSet].Data());
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      
      h_tracks_mean_pt_reso[iSet][iEta]     = (TH1D*)inputFiles[iSet]->Get(Form("histPtResol_%s_%d", addMean.Data() ,iEta));
      h_tracks_sigma_pt_reso[iSet][iEta]    = (TH1D*)inputFiles[iSet]->Get(Form("histPtResol_%s_%d", addSigma.Data() ,iEta));
      
      h_tracks_mean_pt_resoEta[iSet][iEta]  = (TH1D*)inputFiles[iSet]->Get(Form("histEtaResol_%s_%d", addMean.Data() ,iEta));
      h_tracks_sigma_pt_resoEta[iSet][iEta] = (TH1D*)inputFiles[iSet]->Get(Form("histEtaResol_%s_%d", addSigma.Data() ,iEta));

      h_tracks_mean_pt_resoPhi[iSet][iEta]  = (TH1D*)inputFiles[iSet]->Get(Form("histPhiResol_%s_%d", addMean.Data() ,iEta));
      h_tracks_sigma_pt_resoPhi[iSet][iEta] = (TH1D*)inputFiles[iSet]->Get(Form("histPhiResol_%s_%d", addSigma.Data() ,iEta));

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
    DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.02, 0.105);
    // cReso->SetLogz();

    TH2F* histoDummyPtResMean   = new TH2F("histoDummyPtResMean","histoDummyPtResMean",1000,0, 10,1000,-0.1, 0.1);
    SetStyleHistoTH2ForGraphs(histoDummyPtResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPtResMean->GetXaxis()->SetNoExponent();
    histoDummyPtResMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPtResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoDummyPtResMean->Draw();
    TLegend* legendPtResM  = GetAndSetLegend2(0.14, 0.94-(nSets*textSizeLabelsRel), 0.35, 0.94,1.1*textSizeLabelsPixel, 1, "", 43, 0.2);
    DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
    
    for(Int_t iSet=0; iSet<nSets;iSet++){
      DrawGammaSetMarker(h_tracks_mean_pt_reso[iSet][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      h_tracks_mean_pt_reso[iSet][iEta]->Draw("same,p");
      legendPtResM->AddEntry(h_tracks_mean_pt_reso[iSet][iEta],labels[iSet].Data(),"p");
    }
    legendPtResM->Draw();
    drawLatexAdd(Form("#pi^{-} in %1.1f<#eta<%1.1f",etaMin,etaMax),0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PtResolution_Mean_pT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    TH2F* histoDummyPtResSigma   = new TH2F("histoDummyPtResSigma","histoDummyPtResSigma",1000,0, 10,1000,-0.0, maxPtSigma);
    SetStyleHistoTH2ForGraphs(histoDummyPtResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPtResSigma->GetXaxis()->SetNoExponent();
    histoDummyPtResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPtResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoDummyPtResSigma->Draw();
    
    for(Int_t iSet=0; iSet<nSets;iSet++){
      DrawGammaSetMarker(h_tracks_sigma_pt_reso[iSet][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      h_tracks_sigma_pt_reso[iSet][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    drawLatexAdd(Form("#pi^{-} in %1.1f<#eta<%1.1f",etaMin,etaMax),0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PtResolution_Sigma_pT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  }
//   DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.045, 0.105);
//   TH2F* histoDummyEtaResMean   = new TH2F("histoDummyEtaResMean","histoDummyEtaResMean",1000,0, 10,1000,-0.001, 0.001);
//   SetStyleHistoTH2ForGraphs(histoDummyEtaResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#eta^{rec} - #eta^{MC}) / #eta^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
//   histoDummyEtaResMean->GetXaxis()->SetNoExponent();
//   histoDummyEtaResMean->GetYaxis()->SetNdivisions(510,kTRUE);
//   histoDummyEtaResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histoDummyEtaResMean->Draw();
//   TLegend* legendEtaResM  = GetAndSetLegend2(0.14, 0.91-(nActiveEta/2*textSizeLabelsRel), 0.55, 0.91,1.1*textSizeLabelsPixel, 2, "", 43, 0.2);
//   DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
//   
//   for(Int_t iEta=0; iEta<nEta+1;iEta++){
//     if (!enablePlot[iEta]) continue;
//     DrawGammaSetMarker(h_tracks_mean_pt_resoEta[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
//     h_tracks_mean_pt_resoEta[iEta]->Draw("same,p");
//     if (iEta == nEta)
//       legendEtaResM->AddEntry(h_tracks_mean_pt_resoEta[iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
//     else 
//       legendEtaResM->AddEntry(h_tracks_mean_pt_resoEta[iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
//   }
//   legendEtaResM->Draw();
//   drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
//   cReso->Print(Form("%s/EtaResolution_Mean_pT.%s", outputDir.Data(), suffix.Data()));
// 
//   
//   TH2F* histoDummyEtaResSigma   = new TH2F("histoDummyEtaResSigma","histoDummyEtaResSigma",1000,0, 10,1000,-0.0, maxEtaSigma);
//   SetStyleHistoTH2ForGraphs(histoDummyEtaResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#eta^{rec} - #eta^{MC}) / #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
//   histoDummyEtaResSigma->GetXaxis()->SetNoExponent();
//   histoDummyEtaResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
//   histoDummyEtaResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histoDummyEtaResSigma->Draw();
//   
//   for(Int_t iEta=0; iEta<nEta+1;iEta++){
//     if (!enablePlot[iEta]) continue;
//     DrawGammaSetMarker(h_tracks_sigma_pt_resoEta[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
//     h_tracks_sigma_pt_resoEta[iEta]->Draw("same,p");
//   }
//   legendEtaResM->Draw();
//   drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
//   cReso->Print(Form("%s/EtaResolution_Sigma_pT.%s", outputDir.Data(), suffix.Data()));
// 
//   TH2F* histoDummyPhiResMean   = new TH2F("histoDummyPhiResMean","histoDummyPhiResMean",1000,0, 10,1000,-0.01, 0.01);
//   SetStyleHistoTH2ForGraphs(histoDummyPhiResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#varphi^{rec} - #varphi^{MC}) / #varphi^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
//   histoDummyPhiResMean->GetXaxis()->SetNoExponent();
//   histoDummyPhiResMean->GetYaxis()->SetNdivisions(510,kTRUE);
//   histoDummyPhiResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histoDummyPhiResMean->Draw();
//   TLegend* legendPhiResM  = GetAndSetLegend2(0.14, 0.91-(nActiveEta/2*textSizeLabelsRel), 0.55, 0.91,1.1*textSizeLabelsPixel, 2, "", 43, 0.2);
//   DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
//   
//   for(Int_t iEta=0; iEta<nEta+1;iEta++){
//     if (!enablePlot[iEta]) continue;
//     DrawGammaSetMarker(h_tracks_mean_pt_resoPhi[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
//     h_tracks_mean_pt_resoPhi[iEta]->Draw("same,p");
//     if (iEta == nEta)
//       legendPhiResM->AddEntry(h_tracks_mean_pt_resoPhi[iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
//     else 
//       legendPhiResM->AddEntry(h_tracks_mean_pt_resoPhi[iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
//   }
//   legendPhiResM->Draw();
//   drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
//   cReso->Print(Form("%s/PhiResolution_Mean_pT.%s", outputDir.Data(), suffix.Data()));
// 
//   TH2F* histoDummyPhiResSigma   = new TH2F("histoDummyPhiResSigma","histoDummyPhiResSigma",1000,0, 10,1000,-0.0, maxPhiSigma);
//   SetStyleHistoTH2ForGraphs(histoDummyPhiResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#varphi^{rec} - #varphi^{MC}) / #varphi^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
//   histoDummyPhiResSigma->GetXaxis()->SetNoExponent();
//   histoDummyPhiResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
//   histoDummyPhiResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histoDummyPhiResSigma->Draw();
//   
//   for(Int_t iEta=0; iEta<nEta+1;iEta++){
//     if (!enablePlot[iEta]) continue;
//     DrawGammaSetMarker(h_tracks_sigma_pt_resoPhi[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
//     h_tracks_sigma_pt_resoPhi[iEta]->Draw("same,p");
//   }
//   legendPhiResM->Draw();
//   drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
//   cReso->Print(Form("%s/PhiResolution_Sigma_pT.%s", outputDir.Data(), suffix.Data()));
  

}
