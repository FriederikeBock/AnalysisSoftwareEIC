// #include <Riostream.h>
// #include "TMath.h"
// #include <stdlib.h>
// #include <fstream>
// #include <math.h>
// #include <TROOT.h>
// #include <TApplication.h>
// #include <TPaveLabel.h>
// #include <TSystem.h>
// #include <TFrame.h>
// #include <TStyle.h>
// #include <TString.h>
// #include "TGaxis.h"
// #include "TFile.h"
// #include "TH1F.h"
// #include "TH1D.h"
// #include "TH2F.h"
// #include "TF1.h"
// #include "TVirtualFitter.h"
// #include "TObject.h"
// #include "TCanvas.h"
// #include "TMultiGraph.h"
// #include "TLegend.h"
// #include "TDatabasePDG.h"
// #include "TMinuit.h"
// #include "TBenchmark.h"
// #include "TRandom.h"
// #include "TLatex.h"
// #include "TASImage.h"
// #include "TPostScript.h"
// #include "TGraphErrors.h"
// #include "TArrow.h"
// #include "TGraphAsymmErrors.h"
// #include "TGaxis.h"
// #include "TMarker.h"
// #include "Math/WrappedTF1.h"
// #include "Math/BrentRootFinder.h"
// #include "TFitResultPtr.h"
// #include "TFitResult.h"
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value,const Float_t min_towers, const bool normalize);
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);


void resolutionHCAL(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "e-p, #sqrt{#it{s}} = 100 GeV";

  TString outputDir 						                  = Form("plotsResolutionHCAL");
  gSystem->Exec("mkdir -p "+outputDir);


  //************************** Read data **************************************************
  const Int_t nEne = 25;
  const static Double_t partE[]   = {1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0,  10., 11., 13., 15., 20., 25., 30., 40., 50., 60., 70., 80., 150.};
  Int_t useE[nEne]                = {  0,   0,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1};
  Int_t useEgamma[nEne]           = {  0,   0,   0,   1,   1,   1,   0,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,   1};
  Double_t fitrange_low[nEne]     = {0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6,0.65, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7};
  Double_t fitrange_hig[nEne]     = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0};
  TH1F* h_pion_fhcal_E[nEne] 	    = {NULL};
  TH1F* h_pion_fhcal_E_allclus[nEne] 	    = {NULL};
  TH1F* h_gamma_fhcal_E[nEne] 	    = {NULL};
  Int_t activeE = 0;
  Int_t activeEgamma = 0;

  Int_t useEBoth[nEne]            = {  0,   0,   0,   0,   0,   1,   0,   1,   1,   0,   1,   0,   1,   0,   0,   1,   1,   0,   1,   1,   0,   0,   0,   0,   0};
  TH1F* h_pion_ecalfhcal_E[nEne] 	    = {NULL};
  Int_t activeEboth = 0;


  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    cout << partE[epart] << endl;
    if(useE[epart]){
      TTree *const t_pion_fhcal =  (TTree *) (new TFile(Form("/home/nschmidt/fun4all_ForwardAna/EICdetector/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_cluster");
      if(t_pion_fhcal){
        activeE++;
        h_pion_fhcal_E[epart] 	= new TH1F(Form("h_pion_fhcal_E_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_allclus[epart] 	= new TH1F(Form("h_pion_fhcal_E_allclus_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_fhcal_E[epart], t_pion_fhcal, partE[epart], 0.3, 2, true);
        fill_histogram(h_pion_fhcal_E_allclus[epart], t_pion_fhcal, partE[epart], 0.3, 1, true);
      }
    }
    if(useEgamma[epart]){

      TTree *const t_gamma_fhcal =  (TTree *) (new TFile(Form("/home/nschmidt/fun4all_ForwardAna/FHCAL_gamma/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_cluster");
      if(t_gamma_fhcal){
        activeEgamma++;
        h_gamma_fhcal_E[epart] 	= new TH1F(Form("h_gamma_fhcal_E_%d",epart), "", 200, 0.0, 2.0);

        fill_histogram(h_gamma_fhcal_E[epart], t_gamma_fhcal, partE[epart], 0.3, 2, true);
      }
    }
    if(useEBoth[epart]){
      TTree *const t_pion_fecalfhcal =  (TTree *) (new TFile(Form("/home/nschmidt/fun4all_ForwardAna/FECALandFHCAL/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_cluster");
      if(t_pion_fecalfhcal){
        activeEboth++;
        h_pion_ecalfhcal_E[epart] 	= new TH1F(Form("h_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_ecalfhcal_E[epart], t_pion_fecalfhcal, partE[epart], 0.3, 2, true);
      }
    }
  }

  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
	split_canvas(cPNG, "cPNG1", activeE);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 1.5,1000,0, 0.070);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_reso[nEne] = {NULL};
  TF1 *fit_reso_NCell[nEne] = {NULL};
  TF1 *fit_reso_gamma[nEne] = {NULL};

  TH1F*  h_pion_resolution_fhcal_E 	= new TH1F("h_pion_resolution_fhcal_E", "", 400, 0.0, 200.0);
  TH1F*  h_pion_resolution_fhcal_1oE 	= new TH1F("h_pion_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  Int_t padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummy->DrawCopy();
    DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
    h_pion_fhcal_E[epart]->Draw("same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    fit_reso[epart] = new TF1(Form("fit_reso_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);

    h_pion_fhcal_E[epart]->Fit(fit_reso[epart],"L0RMEQ");
    fit_reso[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
    fit_reso[epart]->Draw("same");

    drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    h_pion_resolution_fhcal_E->SetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(2)/fit_reso[epart]->GetParameter(2),2)+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2))); //+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)
    h_pion_resolution_fhcal_1oE->SetBinContent(h_pion_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_resolution.%s", outputDir.Data(), suffix.Data()));
  
  TH2F * histEPlotDummyGammaPion                           = new TH2F("histEPlotDummyGammaPion","histEPlotDummyGammaPion",1000,0, 1.5,1000,0, 0.13);
  SetStyleHistoTH2ForGraphs(histEPlotDummyGammaPion, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  histEPlotDummyGammaPion->GetXaxis()->SetNoExponent();
  histEPlotDummyGammaPion->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyGammaPion->GetXaxis()->SetMoreLogLabels(kTRUE);

  split_canvas(cPNG, "cPNG2", activeEgamma);

  TH1F*  h_gamma_resolution_fhcal_E 	= new TH1F("h_gamma_resolution_fhcal_E", "", 400, 0.0, 200.0);
  TH1F*  h_gamma_resolution_fhcal_1oE 	= new TH1F("h_gamma_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart] && !useEgamma[epart]) continue;
    if(!useEgamma[epart]) continue;
    if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyGammaPion->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV}",partE[epart]),0.90,0.22,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("HCAL",0.90,0.16,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart]){
      DrawGammaSetMarker(h_pion_fhcal_E_allclus[epart], 24, 2, kGray+1, kGray+1);
      h_pion_fhcal_E_allclus[epart]->Draw("same");
      DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      h_pion_fhcal_E[epart]->Draw("same");
      DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      fit_reso[epart]->Draw("same");
      // fit_reso[epart]->SetRange(0.6,1.2);
      // DrawGammaSetMarkerTF1(fit_reso[epart], 2, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd("#color[602]{#pi^{-}}",0.25,0.60,2.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    if(useEgamma[epart]){
      DrawGammaSetMarker(h_gamma_fhcal_E[epart], 20, 2, kGreen+2, kGreen+2);
      h_gamma_fhcal_E[epart]->Draw("same");

      fit_reso_gamma[epart] = new TF1(Form("fit_reso_%d",epart), "gaus", 0.75,1.0);
      h_gamma_fhcal_E[epart]->Fit(fit_reso_gamma[epart],"L0RMEQ");
      fit_reso_gamma[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_gamma[epart], 1, 2, kBlack);
      fit_reso_gamma[epart]->Draw("same");

      drawLatexAdd(Form("mean: %1.2f ",fit_reso_gamma[epart]->GetParameter(1)),0.9,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gamma[epart]->GetParameter(2)),0.9,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gamma[epart]->GetParameter(2)/fit_reso_gamma[epart]->GetParameter(1)),0.9,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#color[418]{#gamma}",0.85,0.60,2.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      h_gamma_resolution_fhcal_E->SetBinContent(h_gamma_resolution_fhcal_E->FindBin(partE[epart]),100*fit_reso_gamma[epart]->GetParameter(2)/fit_reso_gamma[epart]->GetParameter(1));
      h_gamma_resolution_fhcal_1oE->SetBinContent(h_gamma_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso_gamma[epart]->GetParameter(2)/fit_reso_gamma[epart]->GetParameter(1));
    }
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_gamma_pion_resolution.%s", outputDir.Data(), suffix.Data()));
  
  split_canvas(cPNG, "cPNG3", activeE);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart] ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyGammaPion->GetYaxis()->SetRangeUser(0.,0.099);
    histEPlotDummyGammaPion->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV}",partE[epart]),0.90,0.22,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{-} on HCAL",0.90,0.16,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart]){
      DrawGammaSetMarker(h_pion_fhcal_E_allclus[epart], 24, 2, kGray+1, kGray+1);
      h_pion_fhcal_E_allclus[epart]->Draw("same");

      fit_reso_NCell[epart] = new TF1(Form("fit_reso_NCell_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);

      h_pion_fhcal_E_allclus[epart]->Fit(fit_reso_NCell[epart],"L0RMEQ");
      fit_reso_NCell[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_NCell[epart], 1, 2, kGray+3);
      fit_reso_NCell[epart]->Draw("same");

      drawLatexAdd(Form("mean: %1.2f ",fit_reso_NCell[epart]->GetParameter(1)),0.9,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso_NCell[epart]->GetParameter(2)),0.9,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_NCell[epart]->GetParameter(2)/fit_reso_NCell[epart]->GetParameter(1)),0.9,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      h_pion_fhcal_E[epart]->Draw("same");
      DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      fit_reso[epart]->Draw("same");
      // fit_reso[epart]->SetRange(0.6,1.2);
      // DrawGammaSetMarkerTF1(fit_reso[epart], 2, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd("#color[923]{all clusters}",0.9,0.60,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#color[602]{#it{N}_{cell} > 1}",0.20,0.60,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_pion_NCell_resolution.%s", outputDir.Data(), suffix.Data()));
  


  TCanvas* cPNGBoth;
  split_canvas(cPNGBoth, "cPNGBoth", activeEboth);

  TH2F * histEPlotDummyBoth                           = new TH2F("histEPlotDummyBoth","histEPlotDummyBoth",1000,0, 1.5,1000,0, 0.045);
  SetStyleHistoTH2ForGraphs(histEPlotDummyBoth, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummyBoth->GetXaxis()->SetNoExponent();
  histEPlotDummyBoth->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyBoth->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_resoBoth[nEne] = {NULL};

  TH1F*  h_pion_resolution_fecalfhcal_E 	= new TH1F("h_pion_resolution_fecalfhcal_E", "", 200, 0.0, 100.0);
  TH1F*  h_pion_resolution_fecalfhcal_1oE 	= new TH1F("h_pion_resolution_fecalfhcal_1oE", "", 500, 0.0, 1.0);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBoth[epart]) continue;
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    DrawGammaSetMarker(h_pion_ecalfhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
    h_pion_ecalfhcal_E[epart]->Draw("same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} on FECAL+FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    fit_resoBoth[epart] = new TF1(Form("fit_resoBoth_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);

    h_pion_ecalfhcal_E[epart]->Fit(fit_resoBoth[epart],"L0RMEQ");
    fit_resoBoth[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_resoBoth[epart], 1, 2, kBlack);
    fit_resoBoth[epart]->Draw("same");

    drawLatexAdd(Form("mean: %1.2f ",fit_resoBoth[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("sigma: %1.2f ",fit_resoBoth[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    h_pion_resolution_fecalfhcal_E->SetBinContent(h_pion_resolution_fecalfhcal_E->FindBin(partE[epart]),100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1));
    // h_pion_resolution_fecalfhcal_E->SetBinError(h_pion_resolution_fecalfhcal_E->FindBin(partE[epart]),100*TMath::Sqrt(TMath::Power(fit_resoBoth[epart]->GetParError(2)/fit_resoBoth[epart]->GetParameter(2),2)+TMath::Power(fit_resoBoth[epart]->GetParError(1)/fit_resoBoth[epart]->GetParameter(1),2)));
    h_pion_resolution_fecalfhcal_1oE->SetBinContent(h_pion_resolution_fecalfhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1));
    // h_pion_resolution_fecalfhcal_1oE->SetBinError(h_pion_resolution_fecalfhcal_1oE->FindBin(partE[epart]),100*TMath::Sqrt(TMath::Power(fit_resoBoth[epart]->GetParError(2)/fit_resoBoth[epart]->GetParameter(2),2)+TMath::Power(fit_resoBoth[epart]->GetParError(1)/fit_resoBoth[epart]->GetParameter(1),2)));
    padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_HCAL_resolution.%s", outputDir.Data(), suffix.Data()));


  TH2F * histEPlotDummyMaterial                           = new TH2F("histEPlotDummyMaterial","histEPlotDummyMaterial",1000,0, 1.5,1000,0, 0.09);
  SetStyleHistoTH2ForGraphs(histEPlotDummyMaterial, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  histEPlotDummyMaterial->GetXaxis()->SetNoExponent();
  histEPlotDummyMaterial->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyMaterial->GetXaxis()->SetMoreLogLabels(kTRUE);

  split_canvas(cPNG, "cPNG4", activeEboth);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart] && !useEBoth[epart]) continue;
    if(!useEBoth[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyMaterial->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV #pi^{-}}",partE[epart]),0.90,0.26,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart]){
      DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      h_pion_fhcal_E[epart]->Draw("same");
      DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      fit_reso[epart]->Draw("same");
      drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    if(useEBoth[epart]){
      DrawGammaSetMarker(h_pion_ecalfhcal_E[epart], 20, 2, kAzure+6, kAzure+6);
      h_pion_ecalfhcal_E[epart]->Draw("same");

      DrawGammaSetMarkerTF1(fit_resoBoth[epart], 1, 2, kBlack);
      fit_resoBoth[epart]->Draw("same");

      drawLatexAdd(Form("mean: %1.2f ",fit_resoBoth[epart]->GetParameter(1)),0.9,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_resoBoth[epart]->GetParameter(2)),0.9,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1)),0.9,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#color[866]{#font[62]{FHCAL + FECAL}}",0.9,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_pion_material_resolution.%s", outputDir.Data(), suffix.Data()));


// RESOLUTION PLOT

  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
  TH2F * histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 89,1000,0, 39);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
    DrawGammaSetMarker(h_pion_resolution_fhcal_E, 46, 3, kBlack, kBlack);

  h_pion_resolution_fhcal_E->Draw("same,p");
  TF1* fit_reso_E = new TF1("fit_reso_E", "[0]+([1]/TMath::Sqrt(x))+([2]/x)", 0,200);
  fit_reso_E->SetParLimits(0,0.,15);
  // fit_reso_E->SetParLimits(0,9.,10.);
  // fit_reso_E->SetParLimits(1,30.,90);
  fit_reso_E->FixParameter(2,0.08);
  h_pion_resolution_fhcal_E->Fit(fit_reso_E,"RMNE");
  fit_reso_E->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_E, 3, 3, kRed+2);
  fit_reso_E->Draw("same");
  drawLatexAdd("#pi^{-} in FHCAL",0.90,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_E->GetParameter(1),fit_reso_E->GetParameter(0)),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/HCAL_resolution_Fitted.%s", outputDir.Data(), suffix.Data()));

  histResoDummy->Draw();

  DrawGammaSetMarker(h_pion_resolution_fhcal_E, 46, 3, kBlue+2, kBlue+2);
  h_pion_resolution_fhcal_E->SetLineStyle(3);
  h_pion_resolution_fhcal_E->SetLineWidth(3);
  h_pion_resolution_fhcal_E->Draw("same,p");
  DrawGammaSetMarkerTF1(fit_reso_E, 3, 3, kBlue+2);
  fit_reso_E->Draw("same");

  DrawGammaSetMarker(h_pion_resolution_fecalfhcal_E, 24, 3, kBlack, kBlack);
  h_pion_resolution_fecalfhcal_E->SetLineStyle(5);
  h_pion_resolution_fecalfhcal_E->SetLineWidth(3);
  h_pion_resolution_fecalfhcal_E->Draw("same,p");
  TF1* fit_resoboth_E = new TF1("fit_resoboth_E", "[0]+[1]/TMath::Sqrt(x)", 0,200);
  fit_resoboth_E->SetParLimits(0,0.,15);
  h_pion_resolution_fecalfhcal_E->Fit(fit_resoboth_E,"QRMNE+");
  fit_resoboth_E->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_resoboth_E, 5, 3, kBlack);
  fit_resoboth_E->Draw("same");

  DrawGammaSetMarker(h_gamma_resolution_fhcal_E, 27, 3, kGreen+2, kGreen+2);
  h_gamma_resolution_fhcal_E->SetLineStyle(4);
  h_gamma_resolution_fhcal_E->SetLineWidth(3);
  h_gamma_resolution_fhcal_E->Draw("same,p");
  TF1* fit_resogamma_E = new TF1("fit_resogamma_E", "[0]+[1]/TMath::Sqrt(x)", 0,200);
  fit_resogamma_E->SetParLimits(0,0.,15);
  h_gamma_resolution_fhcal_E->Fit(fit_resogamma_E,"QRMNE+");
  fit_resogamma_E->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_resogamma_E, 4, 3, kGreen+2);
  fit_resogamma_E->Draw("same");

  TLegend* legendNonLin           = GetAndSetLegend2(0.25, 0.95-(5*textSizeLabelsRel), 0.9, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
  legendNonLin->AddEntry(h_pion_resolution_fhcal_E,Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_E->GetParameter(1),fit_reso_E->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_pion_resolution_fecalfhcal_E,Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_E->GetParameter(1),fit_resoboth_E->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_gamma_resolution_fhcal_E,Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resogamma_E->GetParameter(1),fit_resogamma_E->GetParameter(0)),"pl");
  // legendNonLin->AddEntry((TObject*)0,"","");
  legendNonLin->Draw();

  // drawLatexAdd(Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_E->GetParameter(1),fit_reso_E->GetParameter(0)),0.93,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // drawLatexAdd(Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_E->GetParameter(1),fit_resoboth_E->GetParameter(0)),0.93,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // drawLatexAdd(Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resogamma_E->GetParameter(1),fit_resogamma_E->GetParameter(0)),0.93,0.70,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/FECALFHCAL_resolution_Fitted.%s", outputDir.Data(), suffix.Data()));

// RESOLUTION PLOT 2

  TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 34);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
    DrawGammaSetMarker(h_pion_resolution_fhcal_1oE, 46, 3, kBlack, kBlack);

  h_pion_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_reso_1oE = new TF1("fit_reso_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_1oE->SetParameter(2,0.3);
  // fit_reso_1oE->SetParLimits(2,0.0,10);
  fit_reso_1oE->FixParameter(2,0.08);
  h_pion_resolution_fhcal_1oE->Fit(fit_reso_1oE,"RMNE");
  fit_reso_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_1oE, 3, 3, kRed+2);
  fit_reso_1oE->Draw("same");
  drawLatexAdd("#pi^{-} in FHCAL",0.90,0.25,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),0.90,0.35,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


  cReso2->Print(Form("%s/HCAL_resolution_1oE_Fitted.%s", outputDir.Data(), suffix.Data()));

  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 39);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
  DrawGammaSetMarker(h_pion_resolution_fhcal_1oE, 46, 3, kBlue+2, kBlue+2);
  h_pion_resolution_fhcal_1oE->SetLineStyle(3);
  h_pion_resolution_fhcal_1oE->SetLineWidth(3);
  h_pion_resolution_fhcal_1oE->Draw("same,p");
  fit_reso_1oE = new TF1("fit_reso_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_1oE->SetParLimits(0,0.,15.);
  fit_reso_1oE->FixParameter(2,0.08);
  h_pion_resolution_fhcal_1oE->Fit(fit_reso_1oE,"RMNE");
  fit_reso_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_1oE, 3, 3, kBlue+2);
  fit_reso_1oE->Draw("same");

  DrawGammaSetMarker(h_pion_resolution_fecalfhcal_1oE, 24, 3, kBlack, kBlack);
  h_pion_resolution_fhcal_1oE->SetLineStyle(5);
  h_pion_resolution_fhcal_1oE->SetLineWidth(3);
  h_pion_resolution_fecalfhcal_1oE->Draw("same,p");
  TF1* fit_resoboth_1oE = new TF1("fit_resoboth_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_resoboth_1oE->SetParLimits(0,0.,15.);
  fit_resoboth_1oE->FixParameter(2,0.08);
  h_pion_resolution_fecalfhcal_1oE->Fit(fit_resoboth_1oE,"RMNE");
  fit_resoboth_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_resoboth_1oE, 3, 3, kBlack);
  fit_resoboth_1oE->Draw("same");

  DrawGammaSetMarker(h_gamma_resolution_fhcal_1oE, 27, 3, kGreen+2, kGreen+2);
  h_pion_resolution_fhcal_1oE->SetLineStyle(7);
  h_pion_resolution_fhcal_1oE->SetLineWidth(3);
  h_gamma_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_reso_gamma_1oE = new TF1("fit_reso_gamma_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_gamma_1oE->SetParLimits(0,0.,15.);
  fit_reso_gamma_1oE->FixParameter(2,0.08);
  h_gamma_resolution_fhcal_1oE->Fit(fit_reso_gamma_1oE,"RMNE");
  fit_reso_gamma_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_gamma_1oE, 3, 3, kGreen+2);
  fit_reso_gamma_1oE->Draw("same");

  // drawLatexAdd(Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_1oE->GetParameter(1),fit_resoboth_1oE->GetParameter(0)),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_1oE->GetParameter(1),fit_reso_gamma_1oE->GetParameter(0)),0.15,0.70,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  legendNonLin           = GetAndSetLegend2(0.15, 0.95-(5*textSizeLabelsRel), 0.7, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
  legendNonLin->AddEntry(h_pion_resolution_fecalfhcal_1oE,Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_1oE->GetParameter(1),fit_resoboth_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_pion_resolution_fhcal_1oE,Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_gamma_resolution_fhcal_1oE,Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_1oE->GetParameter(1),fit_reso_gamma_1oE->GetParameter(0)),"pl");
  // legendNonLin->AddEntry((TObject*)0,"","");
  legendNonLin->Draw();
  cReso2->Print(Form("%s/HCAL_resolution_all_1oE_Fitted.%s", outputDir.Data(), suffix.Data()));

}


void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value, const Float_t min_towers, const bool normalize)
{
	Float_t measured_energy;
//      Float_t true_energy;
	t->SetBranchAddress("e", &measured_energy);
//      t->SetBranchAddress("ge", &true_energy);
	Float_t true_eta;
	t->SetBranchAddress("geta", &true_eta);
	Float_t clusterID;
	t->SetBranchAddress("clusterID", &clusterID);
	Float_t nTowers;
	t->SetBranchAddress("ntowers", &nTowers);
	Float_t eventID;
	t->SetBranchAddress("event", &eventID);

	Int_t nentries = Int_t(t->GetEntries());
  Float_t lastEvt = 0;
	for (Int_t i = 0; i < nentries; ++i) {
		if (t->LoadTree(i) < 0)
			break;

		t->GetEntry(i);
    // if(eventID != lastEvt)cout << eventID << "\tE/Etrue: " << measured_energy/true_energy << endl;
    // else cout << "\tE/Etrue: " << measured_energy/true_energy << endl;
		if (((true_eta > -0.5 && true_eta < 0.5)
      || (true_eta > -3 && true_eta < -2)
      || (true_eta > 1.7 && true_eta < 3.0))  //1.7 < 3
    && (measured_energy > min_value && true_energy > 0.1) // above MIP
    // && clusterID==1
    && nTowers>=min_towers // actual cluster with >1 tower
    ){
			h->Fill(measured_energy / true_energy);
    }
    // if(eventID != lastEvt)lastEvt = eventID;

	}
	if (normalize)
		h->Scale(1 / h->GetEntries());

	h->SetXTitle("E_{cluster} / E_{true}");
	h->SetYTitle("entries / #scale[0.5]{#sum} entries      ");
}

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs)
{
  if(numInputs<10){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(3, 3);
  } else if(numInputs<13){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(4, 3);
  } else if(numInputs<17){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(4, 4);
  } else if(numInputs<21){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(5, 4);
  } else {
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*5);
	  cPNG->Divide(5, 5);
  }

}
