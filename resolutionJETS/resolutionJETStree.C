#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

#include <TROOT.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <iostream>

const int njettypes = 4;
const int firstEtaBin = 7;
const int nInputs = 1;

struct plottingStyleData
{
  Color_t color_jets[njettypes] = {kOrange+2, kMagenta+2, kRed+2}; //, kGreen+2, kBlue+2};
  int marker_jets[njettypes] = {21, 20, 22};//, 29, 47};
  int linestyle_jets[njettypes] = {1, 2, 4};//, 8, 9};
  TString str_jet_type[njettypes] = {"track", "calo" , "full", "all"};
  TString str_jet_type_plot[njettypes] = {"Charged Jets (Track)", "Calo Jets (FHCAL+FEMC)" , "Full Jets (Track+FEMC)", "All Jets"};
  TString collisionSystem = "Pythia 6, e+p: 10#times#kern[0.2]{250} GeV^{2}";
  TString format;
};

void plotting(TH1F *scaleData[nInputs][njettypes][16], TString outputDir, plottingStyleData style, float yMin, float yMax, TString xLabel, TString yLabel);

void resolutionJETStree(
    TString suffix            = "png"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plotsTree/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir+"/Slices");
  gSystem->Exec("mkdir -p "+outputDir+"/EtaResolution");
  gSystem->Exec("mkdir -p "+outputDir+"/EtaScale");
  gSystem->Exec("mkdir -p "+outputDir+"/PhiResolution");
  gSystem->Exec("mkdir -p "+outputDir+"/PhiScale");
  gSystem->Exec("mkdir -p "+outputDir+"/JetEnergyResolution");
  gSystem->Exec("mkdir -p "+outputDir+"/JetEnergyScale");
  gSystem->Exec("mkdir -p "+outputDir+"/JetEfficiency");




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

  plottingStyleData style;
  style.format = suffix;

  TH2F*    histo2D_JES_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH2F*    histo2D_JES_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JES_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JES_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JER_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    histo_JER_pT[nInputs][njettypes][nEta+1] = {{{NULL}}};
 
  TH1F*    h_EtaReso_Width_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    h_EtaReso_Mean_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    h_EtaReso_Mean_Eta[nInputs][njettypes] = {{NULL}};
  TH1F*    h_EtaReso_Width_Eta[nInputs][njettypes] = {{NULL}};

  TH1F*    h_PhiReso_Width_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    h_PhiReso_Mean_E[nInputs][njettypes][nEta+1] = {{{NULL}}};
  TH1F*    h_PhiReso_Mean_Eta[nInputs][njettypes] = {{NULL}};
  TH1F*    h_PhiReso_Width_Eta[nInputs][njettypes] = {{NULL}};

  // Jet Efficiency
  TH1F *h_truth_count[nInputs][njettypes] = {NULL};
  TH1F *h_matched_count[nInputs][njettypes] = {NULL};

  // Load required histograms
  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
        h_EtaReso_Mean_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Mean_%s_Eta", style.str_jet_type[ijr].Data()));
        h_EtaReso_Width_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Width_%s_Eta", style.str_jet_type[ijr].Data()));

        h_PhiReso_Mean_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Mean_%s_Eta", style.str_jet_type[ijr].Data()));
        h_PhiReso_Width_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Width_%s_Eta", style.str_jet_type[ijr].Data()));


        h_truth_count[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_truth_count_%s",  style.str_jet_type[ijr].Data()));
        h_matched_count[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_matched_count_%s",  style.str_jet_type[ijr].Data()));

      for (Int_t eT = 0; eT < nEta+1; eT++){
        histo2D_JES_E[iInp][ijr][eT]	= (TH2F*) inputFiles[iInp]->Get(Form("h_jetscale_%s_E_%d", style.str_jet_type[ijr].Data(),eT));
        histo2D_JES_pT[iInp][ijr][eT]	= (TH2F*) inputFiles[iInp]->Get(Form("h_jetscale_%s_pT_%d", style.str_jet_type[ijr].Data(),eT));

        h_EtaReso_Width_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Width_%s_E_%d", style.str_jet_type[ijr].Data(),eT));
        h_EtaReso_Mean_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Mean_%s_E_%d", style.str_jet_type[ijr].Data(),eT));
        h_PhiReso_Width_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Width_%s_E_%d", style.str_jet_type[ijr].Data(),eT));
        h_PhiReso_Mean_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Mean_%s_E_%d", style.str_jet_type[ijr].Data(),eT));

        histo_JES_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JES_%s_E_%d", style.str_jet_type[ijr].Data(),eT));
        histo_JES_pT[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JES_%s_pT_%d", style.str_jet_type[ijr].Data(),eT));
        histo_JER_E[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JER_%s_E_%d", style.str_jet_type[ijr].Data(),eT));
        histo_JER_pT[iInp][ijr][eT]	= (TH1F*) inputFiles[iInp]->Get(Form("h_JER_%s_pT_%d", style.str_jet_type[ijr].Data(),eT));
        if(!histo_JER_pT[iInp][ijr][eT]) std::cout << Form("h_JER_%s_pT_%d", style.str_jet_type[ijr].Data(),eT) << std::endl;
      }
    }
  }


  // Create projections for slices
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

  // Plot scales
  plotting(histo_JES_E, TString(Form("%s/JetEnergyScale/JES", outputDir.Data())), style, -1, 0.3, TString("#it{E}^{jet}"), TString("Mean((#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}))"));  // Plot jet energy scale
  plotting(h_EtaReso_Mean_E, TString(Form("%s/EtaScale/EtaScale", outputDir.Data())), style, -0.4, 0.4, TString("#it{E}^{jet}"), TString("Mean((#eta^{rec} - #eta^{true}) / #eta^{true}))"));  // Plot jet eta scale
  plotting(h_PhiReso_Mean_E, TString(Form("%s/PhiScale/PhiScale", outputDir.Data())), style, -0.4, 0.4, TString("#it{E}^{jet}"), TString("Mean((#phi^{rec} - #phi^{true}) / #phi^{true}))"));  // Plot jet phi scale

  // Plot resolutions
  plotting(histo_JER_E, TString(Form("%s/JetEnergyResolution/JER_E", outputDir.Data())), style, 0, 1, TString("#it{E}^{jet}"), TString("#sigma(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}"));
  plotting(h_EtaReso_Width_E, TString(Form("%s/EtaResolution/EtaReso", outputDir.Data())), style, 0, 1, TString("#it{E}^{jet}"), TString("#sigma(#eta^{rec} - #eta^{true}) / #eta^{true}"));
  plotting(h_PhiReso_Width_E, TString(Form("%s/PhiResolution/PhiReso", outputDir.Data())), style, 0, 1, TString("#it{E}^{jet}"), TString("#sigma(#phi^{rec} - #phi^{true}) / #phi^{true}"));


// 2D PLOT
//   TCanvas* cSingleSlice = new TCanvas("cSingleSlice","",0,0,1100,1000);
//   DrawGammaCanvasSettings( cSingleSlice, 0.11, 0.01, 0.002, 0.105);
//   // cSingleSlice->SetLogz();

//   TH2F* histoJESSliceDummy   = new TH2F("histoJESSliceDummy","histoJESSliceDummy",1000,-0.99, 0.6,1000,0., 0.1);
//   SetStyleHistoTH2ForGraphs(histoJESSliceDummy, "(#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}","norm. counts ", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
//   histoJESSliceDummy->GetYaxis()->SetNoExponent();
//   histoJESSliceDummy->GetYaxis()->SetNdivisions(505,kTRUE);
//   // histoJESSliceDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

//   if(1){
//     int iInp = 0;
//     int ijr = 1;
//     for (Int_t iEbin=1; iEbin < histo2D_JES_E[iInp][ijr][0]->GetNbinsX(); iEbin++){
//       histoJESSliceDummy->Draw();
//       // DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
//       TLegend* legendJES3  = GetAndSetLegend2(0.35, 0.80-(5*textSizeLabelsRel), 0.6, 0.80-(1*textSizeLabelsRel),1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
//       for(Int_t eT=10; eT<14;eT++){
//         h_projectionPlot[iInp][ijr][eT][iEbin]->Sumw2();
//         h_projectionPlot[iInp][ijr][eT][iEbin]->Rebin(2);
//         DrawGammaSetMarker( h_projectionPlot[iInp][ijr][eT][iEbin], markerStyleEta[eT], 1.5*markerSizeEta[eT], colorEta[eT], colorEta[eT]);
//         h_projectionPlot[iInp][ijr][eT][iEbin]->SetLineWidth(4);
//         h_projectionPlot[iInp][ijr][eT][iEbin]->Draw("same,hist");
//         legendJES3->AddEntry( h_projectionPlot[iInp][ijr][eT][iEbin],Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[eT],partEta[eT+1]),"l");
//       }
//       legendJES3->Draw();
//       drawLatexAdd(collisionSystem.Data(),0.35,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//       drawLatexAdd(Form("anti-#it{k}_{T}, #it{R} = 0.5, FHCAL+FEMC jets"),0.35,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//       drawLatexAdd(Form("%1.1f < #it{E}_{jet} < %1.1f GeV/#it{c}",(200./40)*iEbin,(200./40)*(iEbin+1)),0.35,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//       cSingleSlice->Print(Form("%s/Slices/JES_Slice_Plot_EtaBins%d.%s", outputDir.Data(), iEbin, suffix.Data()));
//     }
//   }


//   TCanvas *canvasEfficiency = new TCanvas("canvasEfficiency","",200,10,1000,900);
//   DrawGammaCanvasSettings(canvasEfficiency,  0.11, 0.01, 0.002, 0.105);
//   for (int i = 0; i < njettypes; i++) {
//     TH1F *efficiency = new TH1F(Form("efficiency_%s", style.str_jet_type[i].Data()), "", 40, 0, 200);
//     SetStyleHistoTH1ForGraphs(efficiency, "#it{E}^{jet}","Jet Efficiency",  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
//     for (int j = 0; j < efficiency->GetSize(); j++) {
//       efficiency->SetBinContent(j, h_matched_count[0][i]->GetBinContent(j) / h_truth_count[0][i]->GetBinContent(j));
//     }
//     DrawGammaSetMarker(efficiency, style.marker_jets[i], 2, style.color_jets[i], style.color_jets[i]);
//     efficiency->DrawCopy();
//     efficiency->Draw("hist");
//     drawLatexAdd(collisionSystem.Data(),0.35,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//     drawLatexAdd(Form("anti-#it{k}_{T}, #it{R} = 0.5, %s jets", style.str_jet_type_plot[i].Data()),0.35,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//     canvasEfficiency->SaveAs(Form("%s/JetEfficiency/JetEfficiency_%s.%s", outputDir.Data(), style.str_jet_type_plot[i].Data(), suffix.Data()));
//     delete efficiency;
//   }

}

void plotting(TH1F *scaleData[nInputs][njettypes][16], TString outputFormat, plottingStyleData style, float yMin, float yMax, TString xLabel, TString yLabel) {
  // SETUP
  TCanvas *canvasJES = new TCanvas("canvasJES", "", 200, 10, 1000, 900);
  DrawGammaCanvasSettings( canvasJES, 0.1, 0.01, 0.01, 0.11);
  Double_t textSizeSinglePad = 0.05;
  Double_t textSizeLabelsPixel = 35;
  Double_t textSizeLabelsRel = 58.0 / 1300;
  TH2F * scaleHist = new TH2F("scaleHist", "scaleHist", 1000, 0, 109, 1000, yMin, yMax);
  SetStyleHistoTH2ForGraphs(scaleHist, xLabel, yLabel, 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  scaleHist->GetXaxis()->SetNoExponent();
  scaleHist->GetYaxis()->SetNdivisions(505,true);
  scaleHist->GetXaxis()->SetMoreLogLabels(true);
  DrawGammaLines(0, 109, 0., 0., 1, kGray+2, 7);
  
  // ITERATE OVER JET TYPES
  for(int ijr=0; ijr<njettypes; ijr++){
    scaleHist->DrawCopy();
    TLegend *scaleLegend = GetAndSetLegend2(0.7, 0.95-((nEta-10)*textSizeLabelsRel), 1.1, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
    // ITERATE OVER ETA RANGES
    for (int eta = firstEtaBin; eta < nEta + 1; eta++) {
      DrawGammaSetMarker(scaleData[0][ijr][eta], markerStyleEta[eta], markerSizeEta[eta], colorEta[eta], colorEta[eta]);

      scaleData[0][ijr][eta]->Draw("same,p");

      if(eta < nEta)
        scaleLegend->AddEntry(scaleData[0][ijr][eta],Form("%1.1f < #it{#eta} < %1.1f",partEta[eta], partEta[eta + 1]),"pl");
      else
        scaleLegend->AddEntry(scaleData[0][ijr][eta],Form("%1.1f < #it{#eta} < %1.1f",1.5, 3.5),"pl");
    }
    // DRAWING AND SAVING
    scaleLegend->Draw();
    drawLatexAdd(Form("%s",style.collisionSystem.Data()), 0.16, 0.90, textSizeLabelsRel, false, false, false);
    drawLatexAdd(Form("#it{p}_{T}^{hard} #geq 5 GeV/#it{c}"), 0.16, 0.85, textSizeLabelsRel, false, false, false);
    drawLatexAdd(Form("anti-k_{T}, #it{R}#kern[0.2]{=}#kern[0.1]{0.5}"), 0.16, 0.80, textSizeLabelsRel, false, false, false);
    scaleHist->Draw("same,axis");
    canvasJES->Print(Form("%s_%s.%s", outputFormat.Data(), style.str_jet_type[ijr].Data(), style.format.Data())); 
  }



}