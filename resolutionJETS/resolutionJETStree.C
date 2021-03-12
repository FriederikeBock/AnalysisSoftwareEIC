#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void resolutionJETStree(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "Pythia 6, e+p: 10#times#kern[0.2]{250} GeV^{2}";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plotsTree/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);


  TFile* inputFile = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/TTLMBPTH5_2021_03_10/jetresolutionhistosSave/output_JRH.root");
  TFile* inputFile_noTTL = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/MBPTH5_2021_03_10/jetresolutionhistosSave/output_JRH.root");
  // TFile* inputFile = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/2021_03_09_clusters_calibrated/jetresolutionhistosSave/output_JRH.root");
  // TFile* inputFile_noTTL = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/2021_03_09_noTTL/jetresolutionhistosSave/output_JRH.root");
  const int njettypes = 5;
  TString str_jet_type[njettypes] = {"track", "full", "hcal", "calo" , "all"};
  TString str_jet_type_plot[njettypes] = {"Track", "Track+FEMC", "FHCal", "FHCal+FEMC" , "Track+FHCal+FEMC (no PF yet)"};
  TH2F*    histo2D_JES_E[njettypes] = {NULL};
  TH2F*    histo2D_JES_pT[njettypes] = {NULL};
  TH1F*    histo_JES_E[njettypes] = {NULL};
  TH1F*    histo_JES_pT[njettypes] = {NULL};
  TH1F*    histo_JER_E[njettypes] = {NULL};
  TH1F*    histo_JER_pT[njettypes] = {NULL};
  TH2F*    histo2D_JES_noTTL_E[njettypes] = {NULL};
  TH2F*    histo2D_JES_noTTL_pT[njettypes] = {NULL};
  TH1F*    histo_JES_noTTL_E[njettypes] = {NULL};
  TH1F*    histo_JES_noTTL_pT[njettypes] = {NULL};
  TH1F*    histo_JER_noTTL_E[njettypes] = {NULL};
  TH1F*    histo_JER_noTTL_pT[njettypes] = {NULL};
  for(int ijr=0;ijr<njettypes;ijr++){
    histo2D_JES_E[ijr]	= (TH2F*) inputFile->Get(Form("h_jetscale_%s_E",str_jet_type[ijr].Data()));
    histo2D_JES_pT[ijr]	= (TH2F*) inputFile->Get(Form("h_jetscale_%s_pT",str_jet_type[ijr].Data()));

    histo_JES_E[ijr]	= (TH1F*) inputFile->Get(Form("h_JES_%s_E",str_jet_type[ijr].Data()));
    histo_JES_pT[ijr]	= (TH1F*) inputFile->Get(Form("h_JES_%s_pT",str_jet_type[ijr].Data()));
    histo_JER_E[ijr]	= (TH1F*) inputFile->Get(Form("h_JER_%s_E",str_jet_type[ijr].Data()));
    histo_JER_pT[ijr]	= (TH1F*) inputFile->Get(Form("h_JER_%s_pT",str_jet_type[ijr].Data()));

    histo2D_JES_noTTL_E[ijr]	= (TH2F*) inputFile_noTTL->Get(Form("h_jetscale_%s_E",str_jet_type[ijr].Data()));
    histo2D_JES_noTTL_pT[ijr]	= (TH2F*) inputFile_noTTL->Get(Form("h_jetscale_%s_pT",str_jet_type[ijr].Data()));

    histo_JES_noTTL_E[ijr]	= (TH1F*) inputFile_noTTL->Get(Form("h_JES_%s_E",str_jet_type[ijr].Data()));
    histo_JES_noTTL_pT[ijr]	= (TH1F*) inputFile_noTTL->Get(Form("h_JES_%s_pT",str_jet_type[ijr].Data()));
    histo_JER_noTTL_E[ijr]	= (TH1F*) inputFile_noTTL->Get(Form("h_JER_%s_E",str_jet_type[ijr].Data()));
    histo_JER_noTTL_pT[ijr]	= (TH1F*) inputFile_noTTL->Get(Form("h_JER_%s_pT",str_jet_type[ijr].Data()));
  }

  Color_t color_jets[njettypes]  = {kOrange+2, kMagenta+2, kRed+2, kGreen+2, kBlue+2};
  int marker_jets[njettypes]     = {21, 20, 22, 29, 47};
  int linestyle_jets[njettypes]  = {1, 2, 4, 8, 9};

  TCanvas *canvasJES                            = new TCanvas("canvasJES","",200,10,1000,900);
  DrawGammaCanvasSettings( canvasJES, 0.1, 0.01, 0.01, 0.11);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histJESPlotDummy_E                           = new TH2F("histJESPlotDummy_E","histJESPlotDummy_E",1000,0, 199,1000,-1.3, 0.89);
  SetStyleHistoTH2ForGraphs(histJESPlotDummy_E, "#it{E}^{jet}","(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histJESPlotDummy_E->GetXaxis()->SetNoExponent();
  histJESPlotDummy_E->GetYaxis()->SetNdivisions(505,kTRUE);
  histJESPlotDummy_E->GetXaxis()->SetMoreLogLabels(kTRUE);

  histJESPlotDummy_E->DrawCopy();
  TLegend* legendJES_E           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  DrawGammaLines(0, 199, 0., 0., 1, kGray+2, 7);
  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(histo_JES_E[ijr], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    histo_JES_E[ijr]->Draw("same,p");

    legendJES_E->AddEntry(histo_JES_E[ijr],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }

  legendJES_E->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histJESPlotDummy_E->Draw("same,axis");

  canvasJES->Print(Form("%s/JES_E_all.%s", outputDir.Data(), suffix.Data()));

  histJESPlotDummy_E->DrawCopy();
  legendJES_E           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  DrawGammaLines(0, 199, 0., 0., 1, kGray+2, 7);
    DrawGammaSetMarker(histo_JES_E[0], marker_jets[0], 2, color_jets[0], color_jets[0]);
    histo_JES_E[0]->Draw("same,p");
    legendJES_E->AddEntry(histo_JES_E[0],"LBL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JES_noTTL_E[0], marker_jets[1], 2, color_jets[1], color_jets[1]);
    histo_JES_noTTL_E[0]->Draw("same,p");
    legendJES_E->AddEntry(histo_JES_noTTL_E[0],"LBL+TTL, #it{B}=1.5 T","pl");

    TH1F* dummJESClone3T_E = (TH1F*)histo_JES_E[0]->Clone("dummJESClone3T_E");
    dummJESClone3T_E->Scale(0.6);
    DrawGammaSetMarker(dummJESClone3T_E, 25, 2, kRed+2, kRed+2);
    dummJESClone3T_E->Draw("same,p");
    legendJES_E->AddEntry(dummJESClone3T_E,"LBL, #it{B}=3.0 T (dummy data)","pl");
    TH1F* dummJESClone3TTTL_E = (TH1F*)histo_JES_noTTL_E[0]->Clone("dummJESClone3TTTL_E");
    dummJESClone3TTTL_E->Scale(0.6);
    DrawGammaSetMarker(dummJESClone3TTTL_E, 24, 2, kBlue+2, kBlue+2);
    dummJESClone3TTTL_E->Draw("same,p");
    legendJES_E->AddEntry(dummJESClone3TTTL_E,"LBL+TTL, #it{B}=3.0 T (dummy data)","pl");


  legendJES_E->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histJESPlotDummy_E->Draw("same,axis");

  canvasJES->Print(Form("%s/JES_E_field_TTL.%s", outputDir.Data(), suffix.Data()));




  TH2F * histJESPlotDummy_pT                           = new TH2F("histJESPlotDummy_pT","histJESPlotDummy_pT",1000,0, 26,1000,-1.3, 0.89);
  SetStyleHistoTH2ForGraphs(histJESPlotDummy_pT, "#it{p}_{T}^{jet}","(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histJESPlotDummy_pT->GetXaxis()->SetNoExponent();
  histJESPlotDummy_pT->GetYaxis()->SetNdivisions(505,kTRUE);
  histJESPlotDummy_pT->GetXaxis()->SetMoreLogLabels(kTRUE);

  histJESPlotDummy_pT->DrawCopy();
  TLegend* legendJES_pT           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  DrawGammaLines(0, 199, 0., 0., 1, kGray+2, 7);
  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(histo_JES_pT[ijr], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    histo_JES_pT[ijr]->Draw("same,p");

    legendJES_pT->AddEntry(histo_JES_pT[ijr],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }

  legendJES_pT->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histJESPlotDummy_pT->Draw("same,axis");

  canvasJES->Print(Form("%s/JES_pT_all.%s", outputDir.Data(), suffix.Data()));

  histJESPlotDummy_pT->DrawCopy();
  legendJES_pT           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  DrawGammaLines(0, 199, 0., 0., 1, kGray+2, 7);
    DrawGammaSetMarker(histo_JES_pT[0], marker_jets[0], 2, color_jets[0], color_jets[0]);
    histo_JES_pT[0]->Draw("same,p");
    legendJES_pT->AddEntry(histo_JES_pT[0],"LBL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JES_noTTL_pT[0], marker_jets[1], 2, color_jets[1], color_jets[1]);
    histo_JES_noTTL_pT[0]->Draw("same,p");
    legendJES_pT->AddEntry(histo_JES_noTTL_pT[0],"LBL+TTL, #it{B}=1.5 T","pl");

    TH1F* dummJESClone3T_pT = (TH1F*)histo_JES_pT[0]->Clone("dummJESClone3T_pT");
    dummJESClone3T_pT->Scale(0.6);
    DrawGammaSetMarker(dummJESClone3T_pT, 25, 2, kRed+2, kRed+2);
    dummJESClone3T_pT->Draw("same,p");
    legendJES_pT->AddEntry(dummJESClone3T_pT,"LBL, #it{B}=3.0 T (dummy data)","pl");
    TH1F* dummJESClone3TTTL_pT = (TH1F*)histo_JES_noTTL_pT[0]->Clone("dummJESClone3TTTL_pT");
    dummJESClone3TTTL_pT->Scale(0.6);
    DrawGammaSetMarker(dummJESClone3TTTL_pT, 24, 2, kBlue+2, kBlue+2);
    dummJESClone3TTTL_pT->Draw("same,p");
    legendJES_pT->AddEntry(dummJESClone3TTTL_pT,"LBL+TTL, #it{B}=3.0 T (dummy data)","pl");


  legendJES_pT->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histJESPlotDummy_pT->Draw("same,axis");

  canvasJES->Print(Form("%s/JES_pT_field_TTL.%s", outputDir.Data(), suffix.Data()));

  TCanvas *canvasJER                            = new TCanvas("canvasJER","",200,10,1000,900);
  DrawGammaCanvasSettings( canvasJER, 0.1, 0.01, 0.01, 0.11);
  textSizeSinglePad                    = 0.05;
  textSizeLabelsPixel                    = 35;
  textSizeLabelsRel                    = 58./1300;
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 199,1000,0, 0.69);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  histEPlotDummy->DrawCopy();
  TLegend* legendHE_pion           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(histo_JER_E[ijr], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    histo_JER_E[ijr]->Draw("same,p");

    legendHE_pion->AddEntry(histo_JER_E[ijr],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }
  legendHE_pion->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.17,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  histEPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JER_E_all.%s", outputDir.Data(), suffix.Data()));



  histEPlotDummy->GetYaxis()->SetRangeUser(0,0.58);
  histEPlotDummy->DrawCopy();
  legendHE_pion           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(4*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

    DrawGammaSetMarker(histo_JER_E[0], marker_jets[0], 2, color_jets[0], color_jets[0]);
    histo_JER_E[0]->Draw("same,p");
    legendHE_pion->AddEntry(histo_JER_E[0],"LBL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JER_noTTL_E[0], marker_jets[1], 2, color_jets[1], color_jets[1]);
    histo_JER_noTTL_E[0]->Draw("same,p");
    legendHE_pion->AddEntry(histo_JER_noTTL_E[0],"LBL+TTL, #it{B}=1.5 T","pl");

    TH1F* dummyClone3T_E = (TH1F*)histo_JER_E[0]->Clone("dummyClone3T_E");
    dummyClone3T_E->Scale(0.6);
    DrawGammaSetMarker(dummyClone3T_E, 25, 2, kRed+2, kRed+2);
    dummyClone3T_E->Draw("same,p");
    legendHE_pion->AddEntry(dummyClone3T_E,"LBL, #it{B}=3.0 T (dummy data)","pl");
    TH1F* dummyClone3TTTL_E = (TH1F*)histo_JER_noTTL_E[0]->Clone("dummyClone3TTTL_E");
    dummyClone3TTTL_E->Scale(0.6);
    DrawGammaSetMarker(dummyClone3TTTL_E, 24, 2, kBlue+2, kBlue+2);
    dummyClone3TTTL_E->Draw("same,p");
    legendHE_pion->AddEntry(dummyClone3TTTL_E,"LBL+TTL, #it{B}=3.0 T (dummy data)","pl");


  legendHE_pion->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.17,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  histEPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JER_B_field_TTL_E.%s", outputDir.Data(), suffix.Data()));

  TH2F * histopTPlotDummy                           = new TH2F("histopTPlotDummy","histopTPlotDummy",1000,0, 26,1000,0, 0.59);
  SetStyleHistoTH2ForGraphs(histopTPlotDummy, "#it{p}_{T}^{jet}","#sigma(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histopTPlotDummy->GetXaxis()->SetNoExponent();
  histopTPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histopTPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  histopTPlotDummy->DrawCopy();
  legendHE_pion           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(histo_JER_pT[ijr], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    histo_JER_pT[ijr]->Draw("same,p");

    legendHE_pion->AddEntry(histo_JER_pT[ijr],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }
  legendHE_pion->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histopTPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JER_pT_all.%s", outputDir.Data(), suffix.Data()));


  histopTPlotDummy->GetYaxis()->SetRangeUser(0,0.37);
  histopTPlotDummy->DrawCopy();
  legendHE_pion           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(4*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

    DrawGammaSetMarker(histo_JER_pT[0], marker_jets[0], 2, color_jets[0], color_jets[0]);
    histo_JER_pT[0]->Draw("same,p");
    legendHE_pion->AddEntry(histo_JER_pT[0],"LBL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JER_noTTL_pT[0], marker_jets[1], 2, color_jets[1], color_jets[1]);
    histo_JER_noTTL_pT[0]->Draw("same,p");
    legendHE_pion->AddEntry(histo_JER_noTTL_pT[0],"LBL+TTL, #it{B}=1.5 T","pl");
    TH1F* dummyClone3T = (TH1F*)histo_JER_pT[0]->Clone("dummyClone3T");
    dummyClone3T->Scale(0.6);
    DrawGammaSetMarker(dummyClone3T, 25, 2, kRed+2, kRed+2);
    dummyClone3T->Draw("same,p");
    legendHE_pion->AddEntry(dummyClone3T,"LBL, #it{B}=3.0 T (dummy data)","pl");
    TH1F* dummyClone3TTTL = (TH1F*)histo_JER_noTTL_pT[0]->Clone("dummyClone3T");
    dummyClone3TTTL->Scale(0.6);
    DrawGammaSetMarker(dummyClone3TTTL, 24, 2, kBlue+2, kBlue+2);
    dummyClone3TTTL->Draw("same,p");
    legendHE_pion->AddEntry(dummyClone3TTTL,"LBL+TTL, #it{B}=3.0 T (dummy data)","pl");


  legendHE_pion->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{1.0}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histopTPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JER_B_field_TTL_pT.%s", outputDir.Data(), suffix.Data()));


}
