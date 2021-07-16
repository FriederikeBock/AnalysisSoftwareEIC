#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

#include <TROOT.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TH2.h>
#include <THStack.h>
#include <TLegend.h>

#include <iostream>

const int njettypes = 5;
const int nInputs = 1;

const int firstEtaBin[njettypes] = {1, 10, 10, 10, 10};
const float min_eta[njettypes]   = {-3.5, 1.5, 1.5, 1.5, 1.5};  // TODO Save this info as metadata...
const float max_eta[njettypes]   = {3.5, 3.5, 3.5, 3.5, 3.5};

struct plottingStyleData
{
  Color_t color_jets[njettypes] = {kOrange+2, kMagenta+2, kBlue+2, kRed+2}; //, kGreen+2, kBlue+2};
  int marker_jets[njettypes] = {21, 20, 20, 22, 29};//, 29, 47};
  int linestyle_jets[njettypes] = {1, 2, 2, 4, 8};//, 8, 9};
  TString str_jet_type[njettypes] = {"track", "calo", "hcal", "nocluster", "emcal"};
  TString str_jet_type_plot[njettypes] = {"Track Jets", "Calo Jets", "HCal Jets", "EMC No Cluster", "EMCal Jets"};
  TString collisionSystem = "Pythia 6, ep 18x100, Q^{2}>100";
  TString jetMatching = "anti-#it{k}_{T}, #it{R} = 0.5";
  TString format;
};

TString *cutString(int jettype) {return new TString(Form("%.1f < #eta < %.1f%s", min_eta[jettype], max_eta[jettype], jettype==0 ? ", pT < 30" : ""));}
void drawInfo(plottingStyleData style, float x, float y, int jettype);
void plotResoOrScale(TH1F *scaleData[nInputs][njettypes][16], TString title, TString outputDir, plottingStyleData style, float yMin, float yMax, TString xLabel, TString yLabel);
void plotSpectra(TH2F *spectra[nInputs][njettypes], plottingStyleData style, TString title, TString outputFormat, TH1F *reco[nInputs][njettypes]=nullptr, TH1F *truth[nInputs][njettypes]=nullptr, float textX=0.22, float textY=0.83);
void plotEfficiency(TH1F *h_matched_count[nInputs][njettypes], TH1F *h_truth_count[nInputs][njettypes], plottingStyleData style, TString outputDir);


void resolutionJETStree(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  gStyle->SetOptStat(0);  //show statistic
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.13);


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
  gSystem->Exec("mkdir -p "+outputDir+"/Spectra");




  TString str_input_type[nInputs] = {"ALL-SI-TTL"};//, "", "", ""};//{"ALLSI-TTL", "ALLSI-noTTL" , "ALLSI-3T", "ALLSI-noTTL-3T","ALLSI-TTL-2xGran"};
  TString str_input_type_plot[nInputs] = {"ALL-SI-TTL"};//, "", "", ""};//{"ALL-SI+TTL (500#mum), BABAR (#it{B}=1.5T)", "ALL-SI, BABAR (#it{B}=1.5T)" , "ALL-SI, BEAST (#it{B}=3.0T)", "ALL-SI+TTL (500#mum), BEAST (#it{B}=3.0T)","ALL-SI+TTL (500#mum), 2x CALO granularity"};
  TFile* inputFiles[nInputs];
  inputFiles[0] = new TFile("/home/tristan/ecce/Singularity/analysis/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/output_JRH.root");

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

  // Spectra
  TH2F *h2D_truth_reco_eta[nInputs][njettypes] = {{NULL}};
  TH1F *h_truth_eta[nInputs][njettypes] = {{NULL}};
  TH1F *h_reco_eta[nInputs][njettypes] = {{NULL}};

  TH2F *h2D_truth_reco_phi[nInputs][njettypes] = {{NULL}};
  TH1F *h_truth_phi[nInputs][njettypes] = {{NULL}};
  TH1F *h_reco_phi[nInputs][njettypes] = {{NULL}};

  TH2F *h2D_truth_reco_E[nInputs][njettypes] = {{NULL}};
  TH1F *h_truth_E[nInputs][njettypes] = {{NULL}};
  TH1F *h_reco_E[nInputs][njettypes] = {{NULL}};

  TH2F *h2D_truth_reco_pT[nInputs][njettypes] = {{NULL}};
  TH1F *h_truth_pT[nInputs][njettypes] = {{NULL}};
  TH1F *h_reco_pT[nInputs][njettypes] = {{NULL}};

  // Jet Efficiency
  TH1F *h_truth_count[nInputs][njettypes] = {NULL};
  TH1F *h_matched_count[nInputs][njettypes] = {NULL};

  TH1 *a = new TH1F();

  // Load required histograms
  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
        h_EtaReso_Mean_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Mean_%s_Eta", style.str_jet_type[ijr].Data()));
        h_EtaReso_Width_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_EtaReso_Width_%s_Eta", style.str_jet_type[ijr].Data()));

        h_PhiReso_Mean_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Mean_%s_Eta", style.str_jet_type[ijr].Data()));
        h_PhiReso_Width_Eta[iInp][ijr]	= (TH1F*) inputFiles[iInp]->Get(Form("h_PhiReso_Width_%s_Eta", style.str_jet_type[ijr].Data()));


        h_truth_count[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_truth_count_%s",  style.str_jet_type[ijr].Data()));
        h_matched_count[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_matched_count_%s",  style.str_jet_type[ijr].Data()));

        h2D_truth_reco_eta[iInp][ijr] = (TH2F*)inputFiles[iInp]->Get(Form("h2D_truth_reco_eta_%s", style.str_jet_type[ijr].Data()));
        h_truth_eta[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_truth_eta_%s", style.str_jet_type[ijr].Data()));
        h_reco_eta[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_reco_eta_%s", style.str_jet_type[ijr].Data()));

        h2D_truth_reco_phi[iInp][ijr] = (TH2F*)inputFiles[iInp]->Get(Form("h2D_truth_reco_phi_%s", style.str_jet_type[ijr].Data()));
        h_truth_phi[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_truth_phi_%s", style.str_jet_type[ijr].Data()));
        h_reco_phi[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_reco_phi_%s", style.str_jet_type[ijr].Data()));

        h2D_truth_reco_E[iInp][ijr] = (TH2F*)inputFiles[iInp]->Get(Form("h2D_truth_reco_E_%s", style.str_jet_type[ijr].Data()));
        h_truth_E[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_truth_E_%s", style.str_jet_type[ijr].Data()));
        h_reco_E[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_reco_E_%s", style.str_jet_type[ijr].Data()));

        h2D_truth_reco_pT[iInp][ijr] = (TH2F*)inputFiles[iInp]->Get(Form("h2D_truth_reco_pT_%s", style.str_jet_type[ijr].Data()));
        h_truth_pT[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_truth_pT_%s", style.str_jet_type[ijr].Data()));
        h_reco_pT[iInp][ijr] = (TH1F*)inputFiles[iInp]->Get(Form("h_reco_pT_%s", style.str_jet_type[ijr].Data()));

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
  TH1D *jesSlices[nInputs][njettypes][nEta+1][40] = {NULL};
  for(int iInp=0;iInp<nInputs;iInp++){
    for(int ijr=0;ijr<njettypes;ijr++){
      for (Int_t eT = 0; eT < nEta+1; eT++){
        for (Int_t i=1; i < histo2D_JES_E[iInp][ijr][eT]->GetNbinsX(); i++){
          jesSlices[iInp][ijr][eT][i] = (TH1D*)histo2D_JES_E[iInp][ijr][eT]->ProjectionY(Form("projectionYdummy%d%d%d%d",iInp,ijr,eT,i), i,i+1,"e");
          jesSlices[iInp][ijr][eT][i]->Scale(1/jesSlices[iInp][ijr][eT][i]->GetEntries());
        }
      }
    }
  }

  // Plot scales
  plotResoOrScale(histo_JES_E, TString("E Scale"), TString(Form("%s/JetEnergyScale/JES", outputDir.Data())), style, -0.6, 0.4, TString("#it{E}^{jet}"), TString("Mean((#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true}))"));  // Plot jet energy scale
  plotResoOrScale(h_EtaReso_Mean_E, TString("#eta Scale"), TString(Form("%s/EtaScale/EtaScale", outputDir.Data())), style, -0.4, 0.4, TString("#it{E}^{jet}"), TString("Mean((#eta^{rec} - #eta^{true}))"));  // Plot jet eta scale
  plotResoOrScale(h_PhiReso_Mean_E, TString("#Phi Scale"), TString(Form("%s/PhiScale/PhiScale", outputDir.Data())), style, -0.4, 0.4, TString("#it{E}^{jet}"), TString("Mean((#Phi^{rec} - #Phi^{true}))"));  // Plot jet phi scale

  // Plot resolutions
  plotResoOrScale(histo_JER_E, TString("E Resolution"), TString(Form("%s/JetEnergyResolution/JER_E", outputDir.Data())), style, 0, 0.6, TString("#it{E}^{jet}"), TString("#sigma((#it{E}^{rec} - #it{E}^{true}) / #it{E}^{true})"));
  plotResoOrScale(h_EtaReso_Width_E, TString("#eta Resolution"), TString(Form("%s/EtaResolution/EtaReso", outputDir.Data())), style, 0, 0.4, TString("#it{E}^{jet}"), TString("#sigma((#eta^{rec} - #eta^{true}))"));
  plotResoOrScale(h_PhiReso_Width_E, TString("#Phi Resolution"), TString(Form("%s/PhiResolution/PhiReso", outputDir.Data())), style, 0, 0.4, TString("#it{E}^{jet}"), TString("#sigma((#Phi^{rec} - #Phi^{true}))"));

  // Plot spectra
  plotSpectra(h2D_truth_reco_eta, style, TString("eta"), TString(Form("%s/Spectra/eta", outputDir.Data())), h_reco_eta, h_truth_eta);
  plotSpectra(h2D_truth_reco_phi, style, TString("phi"), TString(Form("%s/Spectra/phi", outputDir.Data())), h_reco_phi, h_truth_phi, 0.18);
  plotSpectra(h2D_truth_reco_E, style, TString("E"), TString(Form("%s/Spectra/E", outputDir.Data())), h_reco_E, h_truth_E, 0.3);
  plotSpectra(h2D_truth_reco_pT, style, TString("pT"), TString(Form("%s/Spectra/pT", outputDir.Data())), h_reco_pT, h_truth_pT);

  // Plot efficiency
  plotEfficiency(h_matched_count, h_truth_count, style, outputDir);

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
//         jesSlices[iInp][ijr][eT][iEbin]->Sumw2();
//         jesSlices[iInp][ijr][eT][iEbin]->Rebin(2);
//         DrawGammaSetMarker( jesSlices[iInp][ijr][eT][iEbin], markerStyleEta[eT], 1.5*markerSizeEta[eT], colorEta[eT], colorEta[eT]);
//         jesSlices[iInp][ijr][eT][iEbin]->SetLineWidth(4);
//         jesSlices[iInp][ijr][eT][iEbin]->Draw("same,hist");
//         legendJES3->AddEntry( jesSlices[iInp][ijr][eT][iEbin],Form("%1.1f < #it{#eta}_{jet} < %1.1f",partEta[eT],partEta[eT+1]),"l");
//       }
//       legendJES3->Draw();
//       drawLatexAdd(collisionSystem.Data(),0.35,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//       drawLatexAdd(Form("anti-#it{k}_{T}, #it{R} = 0.5, FHCAL+FEMC jets"),0.35,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//       drawLatexAdd(Form("%1.1f < #it{E}_{jet} < %1.1f GeV/#it{c}",(200./40)*iEbin,(200./40)*(iEbin+1)),0.35,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
//       cSingleSlice->Print(Form("%s/Slices/JES_Slice_Plot_EtaBins%d.%s", outputDir.Data(), iEbin, suffix.Data()));
//     }
//   }
}

void plotResoOrScale(TH1F *scaleData[nInputs][njettypes][16], TString title, TString outputFormat, plottingStyleData style, float yMin, float yMax, TString xLabel, TString yLabel) {
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
    if (ijr == 0) {
      scaleHist->GetXaxis()->SetRangeUser(0, 59);
    }
    else {
      scaleHist->GetXaxis()->SetRangeUser(0, 109);
    }
    scaleHist->DrawCopy();
    TLegend *scaleLegend = GetAndSetLegend2(0.7, 0.95-((nEta-firstEtaBin[ijr])*textSizeLabelsRel), 1.05, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
    // ITERATE OVER ETA RANGES
    for (int eta = firstEtaBin[ijr]; eta < nEta; eta++) {
      if (eta == nEta - 1) {
        continue; // Skip 3.5-4 range
      }
      DrawGammaSetMarker(scaleData[0][ijr][eta], markerStyleEta[eta], markerSizeEta[eta], colorEta[eta], colorEta[eta]);

      scaleData[0][ijr][eta]->Draw("same,p");

      if(eta < nEta)
        scaleLegend->AddEntry(scaleData[0][ijr][eta],Form("%1.1f < #it{#eta} < %1.1f",partEta[eta], partEta[eta + 1]),"pl");
      else
        scaleLegend->AddEntry(scaleData[0][ijr][eta],Form("%1.1f < #it{#eta} < %1.1f",1.5, 3.5),"pl");
    }
    // DRAWING AND SAVING
    scaleLegend->Draw();
    drawLatexAdd(title.Data(), 0.16, 0.94, textSizeLabelsRel, kFALSE, kFALSE, kFALSE);
    drawInfo(style, 0.16, 0.895, ijr);
    scaleHist->Draw("same,axis");
    canvasJES->Print(Form("%s_%s.%s", outputFormat.Data(), style.str_jet_type[ijr].Data(), style.format.Data())); 
  }
}

void plotSpectra(TH2F *spectra[nInputs][njettypes], plottingStyleData style, TString title, TString outputFormat, TH1F *reco[nInputs][njettypes]=nullptr, TH1F *truth[nInputs][njettypes]=nullptr, float textX=0.22, float textY=0.83) {
  Double_t textSizeLabelsRel = 58.0 / 1300;
  int canvasWidth = 2000;
  int canvasDivisions = 2;
  if (reco != nullptr && truth != nullptr) {
    canvasWidth = 3000;
    canvasDivisions = 3;
  }
  TCanvas *spectraCanvas = new TCanvas(Form("%s spectra", title.Data()), "", 200, 10, canvasWidth, 900);
  spectraCanvas->Divide(canvasDivisions, 1);
  Double_t textSizeSinglePad = 0.05;
  TH1D *recoSpectra[njettypes];
  TH1D *truthSpectra[njettypes];
  THStack *projectionStack[njettypes];
  TLegend *projectionLegend[njettypes];

  for (int i = 0; i < njettypes; i++) {
    spectraCanvas->cd(1);
    SetStyleHistoTH2ForGraphs(spectra[0][i], Form("truth %s", title.Data()), Form("reco %s", title.Data()), 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91, 2);
    // spectra[0][i]->SetTitle(Form("%s, Matched %s, R=0.5", title.Data(), style.str_jet_type_plot[i].Data()));
    gPad->SetLogz();
    spectra[0][i]->Draw("colz");
    spectra[0][i]->GetXaxis()->SetRangeUser(0, 80);
    spectra[0][i]->GetXaxis()->SetNdivisions();
    spectra[0][i]->GetYaxis()->SetRangeUser(0, 80);
    drawInfo(style, textX, textY, i);

    spectraCanvas->cd(2);
    gPad->SetLogy();
    projectionStack[i] = new THStack();
    projectionLegend[i] = new TLegend(.62, .75, 0.85, 0.9); 
    recoSpectra[i] = spectra[0][i]->ProjectionY(Form("reco_spectra_%s_%s", title.Data(), style.str_jet_type[i].Data()));
    recoSpectra[i]->SetTitle("Reco Matched");
    recoSpectra[i]->SetLineColor(kRed);
    projectionStack[i]->Add(recoSpectra[i]);
    projectionLegend[i]->AddEntry(recoSpectra[i], "Reco");

    truthSpectra[i] = spectra[0][i]->ProjectionX(Form("truth_spectra_%s_%s",title.Data(),  style.str_jet_type[i].Data()));
    truthSpectra[i]->SetTitle("Truth Projection");
    truthSpectra[i]->SetLineColor(kBlue);
    truthSpectra[i]->SetLineWidth(4);
    truthSpectra[i]->SetLineStyle(2);
    projectionStack[i]->Add(truthSpectra[i]);
    projectionLegend[i]->AddEntry(truthSpectra[i], "Truth");

    projectionStack[i]->SetMinimum(0.7);
    projectionStack[i]->SetMaximum(projectionStack[i]->GetMaximum() * 5);
    projectionStack[i]->Draw("nostack");
    SetStyleHistoTHStackForGraphs(projectionStack[i], Form("Jet %s", title.Data()), "Counts",  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);


    drawLatexAdd(Form("Matched %s", style.str_jet_type_plot[i].Data()), 0.165, 0.82, 0.056, kFALSE, kFALSE, kFALSE); 
    // projectionStack[i]->SetTitle(Form("Matched %s", style.str_jet_type_plot[i].Data()));
    projectionLegend[i]->SetTextFont(42);
    projectionLegend[i]->Draw();

    THStack *fullStack = new THStack();
    TLegend *fullLegend = new TLegend(.62, .75, 0.85, 0.9); 
    if (reco != nullptr && truth != nullptr) {
      spectraCanvas->cd(3);
      gPad->SetLogy();
      reco[0][i]->SetTitle("Reco full");
      reco[0][i]->SetLineColor(kRed);
      fullStack->Add(reco[0][i]);
      fullLegend->AddEntry(reco[0][i], "Reco");

      truth[0][i]->SetTitle("Truth full");
      truth[0][i]->SetLineColor(kBlue);
      truth[0][i]->SetLineWidth(4);
      truth[0][i]->SetLineStyle(2);
      fullStack->Add(truth[0][i]);
      fullLegend->AddEntry(truth[0][i], "Truth");

      fullStack->SetMinimum(0.7);
      fullStack->SetMaximum(fullStack->GetMaximum() * 5);
      fullStack->Draw("hist nostack");
      SetStyleHistoTHStackForGraphs(fullStack, Form("Jet %s", title.Data()), "Counts",  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
      drawLatexAdd(Form("All %s", style.str_jet_type_plot[i].Data()), 0.17, 0.82, 0.056, kFALSE, kFALSE, kFALSE); 
      // fullStack->SetTitle(Form("All %s", style.str_jet_type_plot[i].Data()));
      fullLegend->SetTextFont(42);
      fullLegend->Draw();
    }


    spectraCanvas->Print(Form("%s_%s.%s", outputFormat.Data(), style.str_jet_type[i].Data(), style.format.Data()));
  }
}

void plotEfficiency(TH1F *h_matched_count[nInputs][njettypes], TH1F *h_truth_count[nInputs][njettypes], plottingStyleData style, TString outputDir) {
  Double_t textSizeSinglePad = 0.05;
  Double_t textSizeLabelsRel = 58.0 / 1300;
  for (int i = 0; i < njettypes; i++) {
    TCanvas *canvasEfficiency = new TCanvas(Form("canvasEfficiency_%s", style.str_jet_type[i].Data()), "", 200, 10, 1000, 900);
    DrawGammaCanvasSettings(canvasEfficiency,  0.11, 0.01, 0.002, 0.105);
    // TH1F *efficiency = new TH1F(Form("efficiency_%s", style.str_jet_type[i].Data()), "", h_truth_count[0][i]->GetNbinsX() - 1, 0, 200);
    TH1F *efficiency = new TH1F(*(h_matched_count[0][i]));
    efficiency->Divide(h_truth_count[0][i]);
    SetStyleHistoTH1ForGraphs(efficiency, "#it{E}^{jet}","Jet Efficiency",  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
    // for (int j = 1; j <= efficiency->GetNbinsX(); j++) {
    //   efficiency->SetBinContent(j, h_matched_count[0][i]->GetBinContent(j) / h_truth_count[0][i]->GetBinContent(j));
    //   if (i == 3) {
    //     std::cout << h_matched_count[0][i]->GetBinContent(j) << "\t" << h_truth_count[0][i]->GetBinContent(j) << "\t" << efficiency->GetBinContent(j) << "\n";
    //   }
    // }
    DrawGammaSetMarker(efficiency, style.marker_jets[i], 2, style.color_jets[i], style.color_jets[i]);
    efficiency->Draw("hist");
    // h_truth_count[0][i]->Draw("hist");
    drawInfo(style, 0.16, 0.90, i);
    canvasEfficiency->SaveAs(Form("%s/JetEfficiency/JetEfficiency_%s.%s", outputDir.Data(), style.str_jet_type_plot[i].Data(), style.format.Data()));
    delete efficiency;
  }
}

void drawInfo(plottingStyleData style, float x, float y, int jettype) {
  Double_t textSizeSinglePad = 0.05;
  Double_t textSizeLabelsRel = 58.0 / 1300;
  drawLatexAdd(style.collisionSystem.Data(), x, y, textSizeLabelsRel, kFALSE, kFALSE, kFALSE);
  drawLatexAdd(Form("%s, %s", style.jetMatching.Data(), style.str_jet_type_plot[jettype].Data()), x, y - 0.05, textSizeLabelsRel, kFALSE, kFALSE, kFALSE);
  drawLatexAdd(cutString(jettype)->Data(), x, y - 0.10, textSizeLabelsRel, false, false, false);
}