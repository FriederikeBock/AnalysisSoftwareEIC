// #include "../common/plottingheader.h"
// #include "../common/binningheader.h"
#include "plottingheader.h"
#include "binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

#include <TROOT.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <iostream>


void resolutionJETStree(
    TString suffix            = "png"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "Pythia 6, e+p: 10#times#kern[0.2]{250} GeV^{2}";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plotsTree/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir+"/Slices");
  gSystem->Exec("mkdir -p "+outputDir+"/EtaResolution");
  gSystem->Exec("mkdir -p "+outputDir+"/PhiResolution");
  gSystem->Exec("mkdir -p "+outputDir+"/JetEnergyResolution");
  gSystem->Exec("mkdir -p "+outputDir+"/JetEnergyScale");
  gSystem->Exec("mkdir -p "+outputDir+"/JetEfficiency");

  const int firstEtaBin = 7;



  const int nInputs = 1;
  TString str_input_type[nInputs] = {"ALL-SI-TTL"};//, "", "", ""};//{"ALLSI-TTL", "ALLSI-noTTL" , "ALLSI-3T", "ALLSI-noTTL-3T","ALLSI-TTL-2xGran"};
  TString str_input_type_plot[nInputs] = {"ALL-SI-TTL"};//, "", "", ""};//{"ALL-SI+TTL (500#mum), BABAR (#it{B}=1.5T)", "ALL-SI, BABAR (#it{B}=1.5T)" , "ALL-SI, BEAST (#it{B}=3.0T)", "ALL-SI+TTL (500#mum), BEAST (#it{B}=3.0T)","ALL-SI+TTL (500#mum), 2x CALO granularity"};
  TFile* inputFiles[nInputs];
  inputFiles[0] = new TFile("/home/tristan/ecce/Singularity/analysis/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/output_JRH.root");
  // inputFiles[0] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250pTHard5/output_JRH.root");
  // inputFiles[1] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/MODULAR_ALLSILICON-TREXTOUT_e10p250pTHard5/output_JRH.root");
  // inputFiles[2] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/BEAST_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh/output_JRH.root");
  // inputFiles[3] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/BEAST_ALLSILICON-TREXTOUT_e10p250pTHard5/output_JRH.root");
  // inputFiles[4] = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT-HC2xEC2x_e10p250pTHard5/output_JRH.root");

  // const Int_t nEta                = 15;
  // Double_t partEta[nEta+1]        = { -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.2, -0.4, 0.4, 1.2,
  //                                      1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

  const int njettypes = 4;
  TString str_jet_type[njettypes] = {"track", "calo" , "full", "all"};
  TString str_jet_type_plot[njettypes] = {"Charged Jets (Track)", "Calo Jets (FHCAL+FEMC)" , "Full Jets (Track+FEMC)", "All Jets"};

  TH2F*    histo2D_JES_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH2F*    histo2D_JES_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JES_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JES_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JER_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JER_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    h_EtaReso_Width_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    h_EtaReso_Mean_Eta[nInputs][njettypes] = {{NULL}};
  TH1F*    h_EtaReso_Width_Eta[nInputs][njettypes] = {{NULL}};
  TH1F*    h_PhiReso_Width_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    h_PhiReso_Mean_Eta[nInputs][njettypes] = {{NULL}};
  TH1F*    h_PhiReso_Width_Eta[nInputs][njettypes] = {{NULL}};

  // Jet Efficiency
  TH1F *h_truth_count[nInputs][njettypes] = {NULL};
  TH1F *h_matched_count[nInputs][njettypes] = {NULL};

  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
        h_EtaReso_Mean_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Mean_%s_Eta",str_jet_type[ijr].Data()));
        h_EtaReso_Width_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Width_%s_Eta",str_jet_type[ijr].Data()));

        h_PhiReso_Mean_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Mean_%s_Eta",str_jet_type[ijr].Data()));
        h_PhiReso_Width_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Width_%s_Eta",str_jet_type[ijr].Data()));


        h_truth_count[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_truth_count_%s", str_jet_type[ijr].Data()));
        h_matched_count[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_matched_count_%s", str_jet_type[ijr].Data()));

      for (Int_t eT = 0; eT < nEta+1; eT++){
        histo2D_JES_E[iInp][ijr][eT]	= (TH2F*) inputFiles[iInp]->Get(Form("h_jetscale_%s_E_%d",str_jet_type[ijr].Data(),eT));
        histo2D_JES_pT[iInp][ijr][eT]	= (TH2F*) inputFiles[iInp]->Get(Form("h_jetscale_%s_pT_%d",str_jet_type[ijr].Data(),eT));

        h_EtaReso_Width_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Width_%s_E_%d",str_jet_type[ijr].Data(),eT));
        h_PhiReso_Width_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Width_%s_E_%d",str_jet_type[ijr].Data(),eT));

        histo_JES_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JES_%s_E_%d",str_jet_type[ijr].Data(),eT));
        histo_JES_pT[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JES_%s_pT_%d",str_jet_type[ijr].Data(),eT));
        histo_JER_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JER_%s_E_%d",str_jet_type[ijr].Data(),eT));
        histo_JER_pT[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JER_%s_pT_%d",str_jet_type[ijr].Data(),eT));
        if(!histo_JER_pT[iInp][ijr][eT]) std::cout << Form("h_JER_%s_pT_%d",str_jet_type[ijr].Data(),eT) << std::endl;
      }
    }
  }



  TH1D *h_projectionPlot[nInputs][njettypes][nEta+1][40] = {NULL};
  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
      for (Int_t eT = 0; eT < nEta+1; eT++){
        for (Int_t i=1; i < histo2D_JES_E[iInp][ijr][eT]->GetNbinsX(); i++){
          h_projectionPlot[iInp][ijr][eT][i] = (TH1D*)histo2D_JES_E[iInp][ijr][eT]->ProjectionY(Form("projectionYdummy%d%d%d%d",iInp,ijr,eT,i), i,i+1,"e");
          h_projectionPlot[iInp][ijr][eT][i]->Scale(1/h_projectionPlot[iInp][ijr][eT][i]->GetEntries());
        }
      }
    }
  }



  int selectedeta = 12;

  Color_t color_jets[njettypes]  = {kOrange+2, kMagenta+2, kRed+2}; //, kGreen+2, kBlue+2};
  int marker_jets[njettypes]     = {21, 20, 22};//, 29, 47};
  int linestyle_jets[njettypes]  = {1, 2, 4};//, 8, 9};

  TCanvas *canvasJES                            = new TCanvas("canvasJES","",200,10,1000,900);
  DrawGammaCanvasSettings( canvasJES, 0.1, 0.01, 0.01, 0.11);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histJESPlotDummy_E                           = new TH2F("histJESPlotDummy_E","histJESPlotDummy_E",1000,0, 109,1000,-0.49, 0.24);
  SetStyleHistoTH2ForGraphs(histJESPlotDummy_E, "#it{E}^{jet}","(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histJESPlotDummy_E->GetXaxis()->SetNoExponent();
  histJESPlotDummy_E->GetYaxis()->SetNdivisions(505,true);
  histJESPlotDummy_E->GetXaxis()->SetMoreLogLabels(true);



  // {"ALLSI-TTL", "ALLSI-noTTL" , "ALLSI-3T", "ALLSI-noTTL-3T"};
  // {"track", "calo" , "all"}
  // eta = 15
  TLegend* legendJES_E           = GetAndSetLegend2(0.16, 0.15, 0.7, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  histJESPlotDummy_E->DrawCopy();


  DrawGammaLines(0, 109, 0., 0., 1, kGray+2, 7);
    // DrawGammaSetMarker(histo_JES_E[1][0][selectedeta], 24, 2, kGray+1, kGray+1);
    // histo_JES_E[1][0][selectedeta]->Draw("same,p");
    // legendJES_E->AddEntry(histo_JES_E[1][0][selectedeta],"ALL-SI, BABAR (#it{B}=1.5T)","pl");

    // DrawGammaSetMarker(histo_JES_E[3][0][selectedeta], 25, 2, kRed-6, kRed-6);
    // histo_JES_E[3][0][selectedeta]->Draw("same,p");
    // legendJES_E->AddEntry(histo_JES_E[3][0][selectedeta],"ALL-SI, BEAST (#it{B}=3.0T)","pl");

    DrawGammaSetMarker(histo_JES_E[0][0][selectedeta], 20, 2, kBlack, kBlack);
    histo_JES_E[0][0][selectedeta]->Draw("same,p");
    legendJES_E->AddEntry(histo_JES_E[0][0][selectedeta],"ALL-SI+TTL (500#mum), BABAR (#it{B}=1.5T)","pl");

    // DrawGammaSetMarker(histo_JES_E[2][0][selectedeta], 21, 2, kRed+2, kRed+2);
    // histo_JES_E[2][0][selectedeta]->Draw("same,p");
    // legendJES_E->AddEntry(histo_JES_E[2][0][selectedeta],"ALL-SI+TTL (500#mum), BEAST (#it{B}=3.0T)","pl");


  legendJES_E->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.95,0.90,textSizeLabelsRel,false,false,true);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.95,0.85,textSizeLabelsRel,false,false,true);
  drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.95,0.80,textSizeLabelsRel,false,false,true);
  if(selectedeta==15)
  drawLatexAdd(Form("1.5 < #it{#eta}_{jet} < 4.0"),0.95,0.75,textSizeLabelsRel,false,false,true);
  else
  drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.95,0.75,textSizeLabelsRel,false,false,true);

  histJESPlotDummy_E->Draw("same,axis");

  canvasJES->Print(Form("%s/JetEnergyResolution/JES_E_field_TTL.%s", outputDir.Data(), suffix.Data()));


  histJESPlotDummy_E                           = new TH2F("histJESPlotDummy_E","histJESPlotDummy_E",1000,0, 109,1000,-0.29, 0.29);
  SetStyleHistoTH2ForGraphs(histJESPlotDummy_E, "#it{E}^{jet}","(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histJESPlotDummy_E->GetXaxis()->SetNoExponent();
  histJESPlotDummy_E->GetYaxis()->SetNdivisions(505,true);
  histJESPlotDummy_E->GetXaxis()->SetMoreLogLabels(true);
  legendJES_E           = GetAndSetLegend2(0.16, 0.15, 0.7, 0.15+(3*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);
  histJESPlotDummy_E->DrawCopy();
  DrawGammaLines(0, 109, 0., 0., 1, kGray+2, 7);
  // histo_JES_E[nInputs][njettypes][nEta+1]
  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(histo_JES_E[0][ijr][selectedeta], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    histo_JES_E[0][ijr][selectedeta]->Draw("same,p");

    legendJES_E->AddEntry(histo_JES_E[0][ijr][selectedeta],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }

  legendJES_E->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.16,0.90,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.16,0.85,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.16,0.80,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%s",str_input_type_plot[0].Data()),0.16,0.75,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.16,0.70,textSizeLabelsRel,false,false,false);

  histJESPlotDummy_E->Draw("same,axis");

  canvasJES->Print(Form("%s/JetEnergyScale/JES_E_all.%s", outputDir.Data(), suffix.Data()));



  // if(histo_JES_E[1][0][selectedeta] && histo_JES_E[1][0][selectedeta]){
  //   TLegend* legend_JES_Bfield           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  //   histJESPlotDummy_E->DrawCopy();
  //   legend_JES_Bfield           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(4*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  //     DrawGammaSetMarker(histo_JES_E[0][0][selectedeta], marker_jets[1], 2, color_jets[1], color_jets[1]);
  //     histo_JES_E[0][0][selectedeta]->Draw("same,p");
  //     legend_JES_Bfield->AddEntry(histo_JES_E[0][0][selectedeta],"ALL-SI+TTL","pl");

  //     DrawGammaSetMarker(histo_JES_E[4][0][selectedeta], 24, 2, kBlue+2, kBlue+2);
  //     histo_JES_E[4][0][selectedeta]->Draw("same,p");
  //     legend_JES_Bfield->AddEntry(histo_JES_E[4][0][selectedeta],"ALL-SI+TTL, 2x CALO granularity","pl");


  //   legend_JES_Bfield->Draw();
  //   drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,false,false,false);
  //   drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,false,false,false);
  //   drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.17,0.80,textSizeLabelsRel,false,false,false);
  // if(selectedeta==15)
  // drawLatexAdd(Form("1.5 < #it{#eta}_{jet} < 4.0"),0.17,0.75,textSizeLabelsRel,false,false,false);
  // else
  // drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.17,0.75,textSizeLabelsRel,false,false,false);

  //   histJESPlotDummy_E->Draw("same,axis");

  //   canvasJES->Print(Form("%s/JetEnergyScale/JES_Granularity_TTL_E.%s", outputDir.Data(), suffix.Data()));
  // }

  TCanvas *canvasJER                            = new TCanvas("canvasJER","",200,10,1000,900);
  DrawGammaCanvasSettings( canvasJER, 0.1, 0.01, 0.01, 0.11);
  textSizeSinglePad                    = 0.05;
  textSizeLabelsPixel                    = 35;
  textSizeLabelsRel                    = 58./1300;
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.39);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);

  histEPlotDummy->DrawCopy();
  TLegend* legend_JER_Bfield           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(5*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  histEPlotDummy->GetYaxis()->SetRangeUser(0,0.58);
  histEPlotDummy->DrawCopy();
  legend_JER_Bfield           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(4*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);
    DrawGammaSetMarker(histo_JER_E[1][0][selectedeta],24, 2, kGray+1, kGray+1);
    histo_JER_E[1][0][selectedeta]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[1][0][selectedeta],"ALL-SI, BABAR (#it{B}=1.5T)","pl");

    DrawGammaSetMarker(histo_JER_E[3][0][selectedeta], 25, 2, kRed-6, kRed-6);
    histo_JER_E[3][0][selectedeta]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[3][0][selectedeta],"ALL-SI, BEAST (#it{B}=3.0T)","pl");

    DrawGammaSetMarker(histo_JER_E[0][0][selectedeta], 20, 2, kBlack, kBlack);
    histo_JER_E[0][0][selectedeta]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[0][0][selectedeta],"ALL-SI+TTL (500#mum), BABAR (#it{B}=1.5T)","pl");

    DrawGammaSetMarker(histo_JER_E[2][0][selectedeta], 21, 2, kRed+2, kRed+2);
    histo_JER_E[2][0][selectedeta]->Draw("same,p");
    legend_JER_Bfield->AddEntry(histo_JER_E[2][0][selectedeta],"ALL-SI+TTL (500#mum), BEAST (#it{B}=3.0T)","pl");


  legend_JER_Bfield->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("charged-only, anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.17,0.80,textSizeLabelsRel,false,false,false);
  if(selectedeta==15)
  drawLatexAdd(Form("1.5 < #it{#eta}_{jet} < 4.0"),0.17,0.75,textSizeLabelsRel,false,false,false);
  else
  drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.17,0.75,textSizeLabelsRel,false,false,false);

  histEPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JetEnergyResolution/JER_B_field_TTL_E.%s", outputDir.Data(), suffix.Data()));


  if(0){
    int ijr = 1;
    histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.55);
    SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
    histEPlotDummy->GetXaxis()->SetNoExponent();
    histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
    histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);
    if(histo_JER_E[0][1][selectedeta] && histo_JER_E[4][1][selectedeta]){
      TLegend* legend_JER_Bfield           = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+(2*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

      histEPlotDummy->DrawCopy();

        DrawGammaSetMarker(histo_JER_E[0][ijr][selectedeta], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);
        histo_JER_E[0][ijr][selectedeta]->Draw("same,p");
        legend_JER_Bfield->AddEntry(histo_JER_E[0][ijr][selectedeta],"ALL-SI+TTL","pl");

        DrawGammaSetMarker(histo_JER_E[4][ijr][selectedeta], 24, 2, kBlue+2, kBlue+2);
        histo_JER_E[4][ijr][selectedeta]->Draw("same,p");
        legend_JER_Bfield->AddEntry(histo_JER_E[4][ijr][selectedeta],"ALL-SI+TTL, 2x CALO granularity","pl");


      legend_JER_Bfield->Draw();
      drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}, %s",str_jet_type_plot[ijr].Data()),0.17,0.80,textSizeLabelsRel,false,false,false);
    if(selectedeta==15)
    drawLatexAdd(Form("1.5 < #it{#eta}_{jet} < 4.0"),0.17,0.75,textSizeLabelsRel,false,false,false);
    else
    drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.17,0.75,textSizeLabelsRel,false,false,false);

      histEPlotDummy->Draw("same,axis");

      canvasJER->Print(Form("%s/JetEnergyResolution/JER_Granularity_TTL_E.%s", outputDir.Data(), suffix.Data()));
    }
  }
  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.79);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);
  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
      histEPlotDummy->DrawCopy();
      TLegend* legendJER_allEta           = GetAndSetLegend2(0.7, 0.95-((nEta-10)*textSizeLabelsRel), 1.1, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
      for (Int_t eT = firstEtaBin; eT < nEta+1; eT++){
        if(eT==14) continue;
        DrawGammaSetMarker(histo_JER_E[iInp][ijr][eT], markerStyleEta[eT], markerSizeEta[eT], colorEta[eT], colorEta[eT]);
        histo_JER_E[iInp][ijr][eT]->Draw("same,p");
        if(eT<nEta)
          legendJER_allEta->AddEntry(histo_JER_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",partEta[eT], partEta[eT+1]),"pl");
        else
          legendJER_allEta->AddEntry(histo_JER_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",1.5, 3.5),"pl");
      }
      legendJER_allEta->Draw();
      drawLatexAdd(Form("%s",collisionSystem.Data()),0.17,0.90,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.17,0.85,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.17,0.80,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("%s",str_jet_type_plot[ijr].Data()),0.17,0.75,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("%s",str_input_type_plot[iInp].Data()),0.17,0.70,textSizeLabelsRel,false,false,false);
      histEPlotDummy->Draw("same,axis");
      canvasJER->Print(Form("%s/JetEnergyResolution/JER_E_%s_%s.%s", outputDir.Data(), str_input_type[iInp].Data(), str_jet_type[ijr].Data(), suffix.Data()));
    }
  }

  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.59);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma((#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);

  histEPlotDummy->DrawCopy();
  TLegend* legendHE_pion           = GetAndSetLegend2(0.5, 0.75, 0.8, 0.75+(3*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(histo_JER_E[0][ijr][selectedeta], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    histo_JER_E[0][ijr][selectedeta]->Draw("same,p");

    legendHE_pion->AddEntry(histo_JER_E[0][ijr][selectedeta],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }
  legendHE_pion->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.16,0.90,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.16,0.85,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.16,0.80,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%s",str_input_type_plot[0].Data()),0.16,0.75,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.16,0.70,textSizeLabelsRel,false,false,false);

  histEPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/JetEnergyResolution/JER_E_all.%s", outputDir.Data(), suffix.Data()));




// 2D PLOT
  TCanvas* cSingleSlice = new TCanvas("cSingleSlice","",0,0,1100,1000);
  DrawGammaCanvasSettings( cSingleSlice, 0.11, 0.01, 0.002, 0.105);
  // cSingleSlice->SetLogz();

  TH2F* histoJESSliceDummy   = new TH2F("histoJESSliceDummy","histoJESSliceDummy",1000,-0.99, 0.6,1000,0., 0.1);
  SetStyleHistoTH2ForGraphs(histoJESSliceDummy, "(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}","norm. counts ", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  histoJESSliceDummy->GetYaxis()->SetNoExponent();
  histoJESSliceDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  // histoJESSliceDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  if(1){
    int iInp = 0;
    int ijr = 1;
    for (Int_t iEbin=1; iEbin < histo2D_JES_E[iInp][ijr][0]->GetNbinsX(); iEbin++){
      histoJESSliceDummy->Draw();
      // DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
      TLegend* legendJES3  = GetAndSetLegend2(0.35, 0.80-(5*textSizeLabelsRel), 0.6, 0.80-(1*textSizeLabelsRel),1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
      for(Int_t eT=10; eT<14;eT++){
        h_projectionPlot[iInp][ijr][eT][iEbin]->Sumw2();
        h_projectionPlot[iInp][ijr][eT][iEbin]->Rebin(2);
        DrawGammaSetMarker( h_projectionPlot[iInp][ijr][eT][iEbin], markerStyleEta[eT], 1.5*markerSizeEta[eT], colorEta[eT], colorEta[eT]);
        h_projectionPlot[iInp][ijr][eT][iEbin]->SetLineWidth(4);
        h_projectionPlot[iInp][ijr][eT][iEbin]->Draw("same,hist");
        legendJES3->AddEntry( h_projectionPlot[iInp][ijr][eT][iEbin],Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[eT],partEta[eT+1]),"l");
      }
      legendJES3->Draw();
      drawLatexAdd(collisionSystem.Data(),0.35,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("anti-#it{k}_{T}, #it{R} = 0.5, FHCAL+FEMC jets"),0.35,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%1.1f < #it{E}_{jet} < %1.1f GeV/#it{c}",(200./40)*iEbin,(200./40)*(iEbin+1)),0.35,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cSingleSlice->Print(Form("%s/Slices/JES_Slice_Plot_EtaBins%d.%s", outputDir.Data(), iEbin, suffix.Data()));
    }
  }



// 2D PLOT
  TCanvas* cEtaReso = new TCanvas("cEtaReso","",0,0,1100,1000);
  DrawGammaCanvasSettings( cEtaReso, 0.12, 0.01, 0.002, 0.105);
  // cEtaReso->SetLogz();

  TH2F* histoEtaReso   = new TH2F("histoEtaReso","histoEtaReso",1000,0.51, 3.99,1000,0., 0.1);
  SetStyleHistoTH2ForGraphs(histoEtaReso, "#it{#eta}^{jet}","#sigma(#it{#eta}^{rec} - #it{#eta}^{true}) / #it{#eta}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.2);
  histoEtaReso->GetYaxis()->SetNoExponent();
  histoEtaReso->GetYaxis()->SetNdivisions(505,kTRUE);
  // histoEtaReso->GetXaxis()->SetMoreLogLabels(kTRUE);

  if(0){
    int iInp = 0;
    int ijr = 1;
      for(int ijr=0;ijr<njettypes;ijr++){
        histoEtaReso->Draw();
        // DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
        TLegend* legendJES3  = GetAndSetLegend2(0.15, 0.15, 0.6, 0.15+(6*0.9*textSizeLabelsRel),1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
        for(int iInp=0;iInp<nInputs;iInp++){
          DrawGammaSetMarker(h_EtaReso_Width_Eta[iInp][ijr], markerStylePID[iInp], markerSizePID[iInp], colorPID[iInp], colorPID[iInp]);
          h_EtaReso_Width_Eta[iInp][ijr]->Draw("same,p");
          legendJES3->AddEntry(h_EtaReso_Width_Eta[iInp][ijr],Form("%s",str_input_type_plot[iInp].Data()),"pl");
        }
        legendJES3->Draw();
        drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("anti-#it{k}_{T}, #it{R} = 0.5, %s",str_jet_type_plot[ijr].Data()),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        // drawLatexAdd(Form("%1.1f < #it{E}_{jet} < %1.1f GeV/#it{c}",(200./40)*iEbin,(200./40)*(iEbin+1)),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        cEtaReso->Print(Form("%s/EtaResolution/EtaResolutionComparison_%s.%s", outputDir.Data(),str_jet_type[ijr].Data(), suffix.Data()));
      }
  }

  canvasJER                            = new TCanvas("canvasJER","",200,10,1000,900);
  DrawGammaCanvasSettings( canvasJER, 0.1, 0.01, 0.01, 0.11);
  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.15);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{#eta}^{rec} - #it{#eta}^{true}) / #it{#eta}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histEPlotDummy->GetYaxis()->SetNoExponent();
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);
  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
      histEPlotDummy->DrawCopy();
      TLegend* legendJER_allEta           = GetAndSetLegend2(0.7, 0.95-((nEta-10)*textSizeLabelsRel), 1.1, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
      for (Int_t eT = firstEtaBin; eT < nEta+1; eT++){
        if(eT==14) continue;
        DrawGammaSetMarker(h_EtaReso_Width_E[iInp][ijr][eT], markerStyleEta[eT], markerSizeEta[eT], colorEta[eT], colorEta[eT]);
        h_EtaReso_Width_E[iInp][ijr][eT]->Draw("same,p");
        if(eT<nEta)
          legendJER_allEta->AddEntry(h_EtaReso_Width_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",partEta[eT], partEta[eT+1]),"pl");
        else
          legendJER_allEta->AddEntry(h_EtaReso_Width_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",1.5, 3.5),"pl");
      }
      legendJER_allEta->Draw();
      drawLatexAdd(Form("%s",collisionSystem.Data()),0.16,0.90-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.16,0.85-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.16,0.80-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("%s",str_jet_type_plot[ijr].Data()),0.16,0.75-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("%s",str_input_type_plot[iInp].Data()),0.16,0.70-.52,textSizeLabelsRel,false,false,false);
      histEPlotDummy->Draw("same,axis");
      canvasJER->Print(Form("%s/EtaResolution/EtaReso_E_%s_%s.%s", outputDir.Data(), str_input_type[iInp].Data(), str_jet_type[ijr].Data(), suffix.Data()));
    }
  }
  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.0599);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{#eta}^{rec} - #it{#eta}^{true}) / #it{#eta}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histEPlotDummy->GetYaxis()->SetNoExponent();
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);

  histEPlotDummy->DrawCopy();
  legendHE_pion           = GetAndSetLegend2(0.5, 0.15, 0.8, 0.15+(3*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(h_EtaReso_Width_E[0][ijr][selectedeta], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    h_EtaReso_Width_E[0][ijr][selectedeta]->Draw("same,p");

    legendHE_pion->AddEntry(h_EtaReso_Width_E[0][ijr][selectedeta],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }
  legendHE_pion->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.16,0.90-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.16,0.85-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.16,0.80-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%s",str_input_type_plot[0].Data()),0.16,0.75-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.16,0.70-.52,textSizeLabelsRel,false,false,false);

  histEPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/EtaResolution/EtaReso_E_all.%s", outputDir.Data(), suffix.Data()));

  // Phi Resolutions
  TCanvas* cPhiReso = new TCanvas("cPhiReso","",0,0,1100,1000);
  DrawGammaCanvasSettings( cPhiReso, 0.12, 0.01, 0.002, 0.105);
  // cPhiReso->SetLogz();

  TH2F* histoPhiReso   = new TH2F("histoPhiReso","histoPhiReso",1000,0.51, 3.99,1000,0., 0.1);
  SetStyleHistoTH2ForGraphs(histoPhiReso, "#it{#phi}^{jet}","#sigma(#it{#phi}^{rec} - #it{#phi}^{true}) / #it{#phi}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.2);
  histoPhiReso->GetYaxis()->SetNoExponent();
  histoPhiReso->GetYaxis()->SetNdivisions(505,kTRUE);
  // histoPhiReso->GetXaxis()->SetMoreLogLabels(kTRUE);

  if(0){
    int iInp = 0;
    int ijr = 1;
      for(int ijr=0;ijr<njettypes;ijr++){
        histoPhiReso->Draw();
        // DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
        TLegend* legendJES3  = GetAndSetLegend2(0.15, 0.15, 0.6, 0.15+(6*0.9*textSizeLabelsRel),1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
        for(int iInp=0;iInp<nInputs;iInp++){
          DrawGammaSetMarker(h_PhiReso_Width_Eta[iInp][ijr], markerStylePID[iInp], markerSizePID[iInp], colorPID[iInp], colorPID[iInp]);
          h_PhiReso_Width_Eta[iInp][ijr]->Draw("same,p");
          legendJES3->AddEntry(h_PhiReso_Width_Eta[iInp][ijr],Form("%s",str_input_type_plot[iInp].Data()),"pl");
        }
        legendJES3->Draw();
        drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("anti-#it{k}_{T}, #it{R} = 0.5, %s",str_jet_type_plot[ijr].Data()),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        // drawLatexAdd(Form("%1.1f < #it{E}_{jet} < %1.1f GeV/#it{c}",(200./40)*iEbin,(200./40)*(iEbin+1)),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        cPhiReso->Print(Form("%s/PhiResolution/PhiResolutionComparison_%s.%s", outputDir.Data(),str_jet_type[ijr].Data(), suffix.Data()));
      }
  }

  canvasJER                            = new TCanvas("canvasJER","",200,10,1000,900);
  DrawGammaCanvasSettings( canvasJER, 0.14, 0.01, 0.01, 0.11);
  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.25);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{#phi}^{rec} - #it{#phi}^{true}) / #it{#phi}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histEPlotDummy->GetYaxis()->SetNoExponent();
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);
  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
      histEPlotDummy->DrawCopy();
      TLegend* legendJER_allEta           = GetAndSetLegend2(0.7, 0.95-((nEta-10)*textSizeLabelsRel), 1.1, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
      for (Int_t eT = firstEtaBin; eT < nEta+1; eT++){
        if(eT==14) continue;
        DrawGammaSetMarker(h_PhiReso_Width_E[iInp][ijr][eT], markerStyleEta[eT], markerSizeEta[eT], colorEta[eT], colorEta[eT]);
        h_PhiReso_Width_E[iInp][ijr][eT]->Draw("same,p");
        if(eT<nEta)
          legendJER_allEta->AddEntry(h_PhiReso_Width_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",partEta[eT], partEta[eT+1]),"pl");
        else
          legendJER_allEta->AddEntry(h_PhiReso_Width_E[iInp][ijr][eT],Form("%1.1f < #it{#eta} < %1.1f",1.5, 3.5),"pl");
      }
      legendJER_allEta->Draw();
      drawLatexAdd(Form("%s",collisionSystem.Data()),0.16,0.90-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.16,0.85-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.16,0.80-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("%s",str_jet_type_plot[ijr].Data()),0.16,0.75-.52,textSizeLabelsRel,false,false,false);
      drawLatexAdd(Form("%s",str_input_type_plot[iInp].Data()),0.16,0.70-.52,textSizeLabelsRel,false,false,false);
      histEPlotDummy->Draw("same,axis");
      canvasJER->Print(Form("%s/PhiResolution/PhiReso_E_%s_%s.%s", outputDir.Data(), str_input_type[iInp].Data(), str_jet_type[ijr].Data(), suffix.Data()));
    }
  }
  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 109,1000,0, 0.0599);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}^{jet}","#sigma(#it{#eta}^{rec} - #it{#eta}^{true}) / #it{#eta}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91, 1, 0.9);
  histEPlotDummy->GetYaxis()->SetNoExponent();
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,true);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(true);

  histEPlotDummy->DrawCopy();
  legendHE_pion           = GetAndSetLegend2(0.5, 0.15, 0.8, 0.15+(3*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

  for(int ijr=0;ijr<njettypes;ijr++){
    DrawGammaSetMarker(h_PhiReso_Width_E[0][ijr][selectedeta], marker_jets[ijr], 2, color_jets[ijr], color_jets[ijr]);

    h_PhiReso_Width_E[0][ijr][selectedeta]->Draw("same,p");

    legendHE_pion->AddEntry(h_PhiReso_Width_E[0][ijr][selectedeta],Form("%s jets",str_jet_type_plot[ijr].Data()),"pl");
  }
  legendHE_pion->Draw();
  drawLatexAdd(Form("%s",collisionSystem.Data()),0.16,0.90-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"),0.16,0.85-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"),0.16,0.80-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%s",str_input_type_plot[0].Data()),0.16,0.75-.52,textSizeLabelsRel,false,false,false);
  drawLatexAdd(Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[selectedeta],partEta[selectedeta+1]),0.16,0.70-.52,textSizeLabelsRel,false,false,false);

  histEPlotDummy->Draw("same,axis");

  canvasJER->Print(Form("%s/PhiResolution/PhiReso_E_all.%s", outputDir.Data(), suffix.Data()));


  TCanvas *canvasEfficiency = new TCanvas("canvasEfficiency","",200,10,1000,900);
  DrawGammaCanvasSettings(canvasEfficiency,  0.11, 0.01, 0.002, 0.105);
  for (int i = 0; i < njettypes; i++) {
    TH1F *efficiency = new TH1F(Form("efficiency_%s", str_jet_type[i].Data()), "", 40, 0, 200);
    SetStyleHistoTH1ForGraphs(efficiency, "#it{E}^{jet}","Jet Efficiency",  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
    for (int j = 0; j < efficiency->GetSize(); j++) {
      efficiency->SetBinContent(j, h_matched_count[0][i]->GetBinContent(j) / h_truth_count[0][i]->GetBinContent(j));
    }
    DrawGammaSetMarker(efficiency, marker_jets[i], 2, color_jets[i], color_jets[i]);
    efficiency->DrawCopy();
    efficiency->Draw("hist");
    drawLatexAdd(collisionSystem.Data(),0.35,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("anti-#it{k}_{T}, #it{R} = 0.5, %s jets", str_jet_type_plot[i].Data()),0.35,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    canvasEfficiency->SaveAs(Form("%s/JetEfficiency/JetEfficiency_%s.%s", outputDir.Data(), str_jet_type_plot[i].Data(), suffix.Data()));
    delete efficiency;
  }

}

