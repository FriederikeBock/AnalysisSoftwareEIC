#include "../common/plottingheader.h"
#include "../common/binningheader.h"
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

  const int nInputs = 4;
  TString str_input_type[nInputs] = {"ALLSI-TTL", "ALLSI-noTTL" , "ALLSI-3T", "ALLSI-noTTL-3T"};
  TString str_input_type_plot[nInputs] = {"ALL-SI+TTL", "ALL-SI" , "ALL-SI+TTL, BEAST (B=3T)", "ALL-SI, BEAST (B=3T)"};
  TFile* inputFiles[nInputs];
  inputFiles[0] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250pTHard5/output_JRH.root");
  inputFiles[1] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/MODULAR_ALLSILICON-TREXTOUT_e10p250pTHard5/output_JRH.root");
  inputFiles[2] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/BEAST_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh/output_JRH.root");
  inputFiles[3] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/BEAST_ALLSILICON-TREXTOUT_e10p250pTHard5/output_JRH.root");

  // const Int_t nEta                = 15;
  // Double_t partEta[nEta+1]        = { -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.2, -0.4, 0.4, 1.2,
  //                                      1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

  const int njettypes = 3;
  TString str_jet_type[njettypes] = {"track", "calo" , "all"};
  TString str_jet_type_plot[njettypes] = {"Charged (Track)", "Calo (FHCAL+FEMC)" , "PF (Track+Calo)"};

  TH2F*    histo2D_JES_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH2F*    histo2D_JES_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JES_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JES_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JER_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JER_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};

  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
      for (Int_t eT = 0; eT < nEta+1; eT++){
        histo2D_JES_E[iInp][ijr][eT]	= (TH2F*) inputFiles[iInp]->Get(Form("h_jetscale_%s_E_%d",str_jet_type[ijr].Data(),eT));
        histo2D_JES_pT[iInp][ijr][eT]	= (TH2F*) inputFiles[iInp]->Get(Form("h_jetscale_%s_pT_%d",str_jet_type[ijr].Data(),eT));

        histo_JES_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JES_%s_E_%d",str_jet_type[ijr].Data(),eT));
        histo_JES_pT[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JES_%s_pT_%d",str_jet_type[ijr].Data(),eT));
        histo_JER_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JER_%s_E_%d",str_jet_type[ijr].Data(),eT));
        histo_JER_pT[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JER_%s_pT_%d",str_jet_type[ijr].Data(),eT));
        if(!histo_JER_pT[iInp][ijr][eT]) cout << Form("h_JER_%s_pT_%d",str_jet_type[ijr].Data(),eT) << endl;
      }
    }
  }

  Color_t color_jets[njettypes]  = {kOrange+2, kMagenta+2, kRed+2}; //, kGreen+2, kBlue+2};
  int marker_jets[njettypes]     = {21, 20, 22}; //, 29, 47};
  int linestyle_jets[njettypes]  = {1, 2, 4}; //, 8, 9};

  TCanvas *canvasJES                            = new TCanvas("canvasJES","",200,10,1000,900);
  DrawGammaCanvasSettings( canvasJES, 0.1, 0.01, 0.01, 0.11);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histJESPlotDummy_E                           = new TH2F("histJESPlotDummy_E","histJESPlotDummy_E",1000,0, 199,1000,-0.49, 0.49);
  SetStyleHistoTH2ForGraphs(histJESPlotDummy_E, "#it{E}^{jet}","(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histJESPlotDummy_E->GetXaxis()->SetNoExponent();
  histJESPlotDummy_E->GetYaxis()->SetNdivisions(505,kTRUE);
  histJESPlotDummy_E->GetXaxis()->SetMoreLogLabels(kTRUE);

  // {"ALLSI-TTL", "ALLSI-noTTL" , "ALLSI-3T", "ALLSI-noTTL-3T"};
  // {"track", "calo" , "all"}
  // eta = 15
  TLegend* legendJES_E           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  histJESPlotDummy_E->DrawCopy();
  legendJES_E           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  DrawGammaLines(0, 199, 0., 0., 1, kGray+2, 7);
    DrawGammaSetMarker(histo_JES_E[1][0][11], marker_jets[0], 2, color_jets[0], color_jets[0]);
    histo_JES_E[1][0][11]->Draw("same,p");
    legendJES_E->AddEntry(histo_JES_E[1][0][11],"LBL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JES_E[3][0][11], 25, 2, kRed+2, kRed+2);
    histo_JES_E[3][0][11]->Draw("same,p");
    legendJES_E->AddEntry(histo_JES_E[3][0][11],"LBL, #it{B}=3.0 T","pl");

    DrawGammaSetMarker(histo_JES_E[0][0][11], marker_jets[1], 2, color_jets[1], color_jets[1]);
    histo_JES_E[0][0][11]->Draw("same,p");
    legendJES_E->AddEntry(histo_JES_E[0][0][11],"LBL+TTL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JES_E[2][0][11], 24, 2, kBlue+2, kBlue+2);
    histo_JES_E[2][0][11]->Draw("same,p");
    legendJES_E->AddEntry(histo_JES_E[2][0][11],"LBL+TTL, #it{B}=3.0 T","pl");


  legendJES_E->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histJESPlotDummy_E->Draw("same,axis");

  canvasJES->Print(Form("%s/JES_E_field_TTL.%s", outputDir.Data(), suffix.Data()));


  // DrawGammaLines(0, 199, 0., 0., 1, kGray+2, 7);
  // for(int ijr=0;ijr<njettypes;ijr++){
  //   DrawGammaSetMarker(histo_JES_E[ijr], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

  //   histo_JES_E[ijr]->Draw("same,p");

  //   legendJES_E->AddEntry(histo_JES_E[ijr],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  // }

  // legendJES_E->Draw();
  // drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  // histJESPlotDummy_E->Draw("same,axis");

  // canvasJES->Print(Form("%s/JES_E_all.%s", outputDir.Data(), suffix.Data()));



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
  TLegend* legend_JER_Bfield           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  histEPlotDummy->GetYaxis()->SetRangeUser(0,0.58);
  histEPlotDummy->DrawCopy();
  legend_JER_Bfield           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(4*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

    DrawGammaSetMarker(histo_JER_E[1][0][11], marker_jets[0], 2, color_jets[0], color_jets[0]);
    histo_JER_E[1][0][11]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[1][0][11],"LBL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JER_E[3][0][11], 25, 2, kRed+2, kRed+2);
    histo_JER_E[3][0][11]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[3][0][11],"LBL, #it{B}=3.0 T","pl");

    DrawGammaSetMarker(histo_JER_E[0][0][11], marker_jets[1], 2, color_jets[1], color_jets[1]);
    histo_JER_E[0][0][11]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[0][0][11],"LBL+TTL, #it{B}=1.5 T","pl");

    DrawGammaSetMarker(histo_JER_E[2][0][11], 24, 2, kBlue+2, kBlue+2);
    histo_JER_E[2][0][11]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[2][0][11],"LBL+TTL, #it{B}=3.0 T","pl");


  legend_JER_Bfield->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.17,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  histEPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JER_B_field_TTL_E.%s", outputDir.Data(), suffix.Data()));

  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
      histEPlotDummy->DrawCopy();
      TLegend* legendJER_allEta           = GetAndSetLegend2(0.7, 0.95-((nEta-10)*textSizeLabelsRel), 1.1, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
      for (Int_t eT = 10; eT < nEta+1; eT++){
        DrawGammaSetMarker(histo_JER_E[iInp][ijr][eT], markerStyleEta[eT], markerSizeEta[eT], colorEta[eT], colorEta[eT]);
        histo_JER_E[iInp][ijr][eT]->Draw("same,p");
        if(eT<nEta)
          legendJER_allEta->AddEntry(histo_JER_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",partEta[eT], partEta[eT+1]),"pl");
        else
          legendJER_allEta->AddEntry(histo_JER_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",1.5, partEta[eT]),"pl");
      }
      legendJER_allEta->Draw();
      drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.17,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%s",str_jet_type_plot[ijr].Data()),0.17,0.75,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%s",str_input_type_plot[iInp].Data()),0.17,0.70,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      histEPlotDummy->Draw("same,axis");
      canvasJER->Print(Form("%s/JER_E_%s_%s.%s", outputDir.Data(), str_input_type[iInp].Data(), str_jet_type[ijr].Data(), suffix.Data()));
    }
  }

/*
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
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

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
  drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.95,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  histopTPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JER_B_field_TTL_pT.%s", outputDir.Data(), suffix.Data()));

*/
}
