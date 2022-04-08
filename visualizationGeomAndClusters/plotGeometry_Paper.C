#include <algorithm>
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
#include "TColor.h"
#include <TChain.h>

void plotGeometry_Paper(
    TString nameFileECal          = "",
    TString nameFileHCal          = "",
    TString suffix            = "pdf"
    
){
  

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  
  TString outputDir                             = Form("plotsPaperGeo");
  gSystem->Exec("mkdir -p "+outputDir);

  double textSizeSinglePad                    = 0.05;
  double textSizeLabelsPixel                  = 35;
  double textSizeLabelsRel                    = 58./1300;

  TFile* inputFileECal                      = new TFile(nameFileECal.Data());
  if(inputFileECal->IsZombie()){
    cout << "inputFileECal not found!" << endl;
  }
  TH2F* h_PhiEta_ECal      = (TH2F*)inputFileECal->Get("h_EtaPhi");
  TH2F* h_PhiEta_iEta_ECal      = (TH2F*)inputFileECal->Get("h_PhiEta_IEta");
  TH2F* h_PhiEta_iPhi_ECal      = (TH2F*)inputFileECal->Get("h_PhiEta_IPhi");

  TFile* inputFileHCal                      = new TFile(nameFileHCal.Data());
  if(inputFileHCal->IsZombie()){
    cout << "inputFileHCal not found!" << endl;
  }
  TH2F* h_PhiEta_HCal      = (TH2F*)inputFileHCal->Get("h_EtaPhi");
  TH2F* h_PhiEta_iEta_HCal      = (TH2F*)inputFileHCal->Get("h_PhiEta_IEta");
  TH2F* h_PhiEta_iPhi_HCal      = (TH2F*)inputFileHCal->Get("h_PhiEta_IPhi");

  TCanvas* cReso = new TCanvas("cReso","",0,0,1000,950);
  
  DrawGammaCanvasSettings( cReso, 0.08, 0.02, 0.01, 0.08);
    SetStyleHistoTH2ForGraphs(h_PhiEta_ECal, "#eta", "#varphi", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.75, 0.73, 510,510);

    h_PhiEta_ECal->SetMaximum(20);
    h_PhiEta_ECal->SetLineColor(kBlack);
    h_PhiEta_ECal->SetFillColor(kBlack);
    h_PhiEta_ECal->Draw("box");
    h_PhiEta_iEta_ECal->SetLineColor(kBlue+2);
    h_PhiEta_iEta_ECal->SetFillColor(kBlue+2);
    h_PhiEta_iEta_ECal->Draw("box,same");
    h_PhiEta_iPhi_ECal->SetLineColor(kRed+2);
    h_PhiEta_iPhi_ECal->SetFillColor(kRed+2);
    h_PhiEta_iPhi_ECal->Draw("box,same");


    TLegend* legendPhiEta    = GetAndSetLegend2(0.12, 0.12, 0.4, 0.12+(3*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendPhiEta->AddEntry(h_PhiEta_ECal, "all towers","plf");
    legendPhiEta->AddEntry(h_PhiEta_iEta_ECal, "towers iEta = 1, 30, 70","fl");
    legendPhiEta->AddEntry(h_PhiEta_iPhi_ECal, "towers iPhi = 1, 50, 100","fl");
    legendPhiEta->Draw();
    
    drawLatexAdd(Form("ECCE"),0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
  cReso->Print(Form("%s/PhiEtaPaper.%s", outputDir.Data(), suffix.Data()));




  Double_t arrayBoundariesX1_4[2];
  Double_t arrayBoundariesY1_4[3];
  Double_t relativeMarginsX[3];
  Double_t relativeMarginsY[3];
  textSizeLabelsPixel             = 50;
  ReturnCorrectValuesForCanvasScaling(1350,1450, 1, 2,0.06, 0.006, 0.005,0.075,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

  TCanvas* canvasGeoPaper     = new TCanvas("canvasGeoPaper","",0,0,1350,1450);  // gives the page size
  DrawGammaCanvasSettings( canvasGeoPaper,  0.13, 0.02, 0.03, 0.06);

  TPad* padGeoECal               = new TPad("padGeoECal", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
  DrawGammaPadSettings( padGeoECal, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
  padGeoECal->Draw();

  TPad* padGeoHCal                = new TPad("padGeoHCal", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
  DrawGammaPadSettings( padGeoHCal, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
  padGeoHCal->Draw();

  padGeoECal->cd();
  // padGeoECal->SetLogx();

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

  TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 800, -4, 4, 400, -TMath::Pi(), TMath::Pi());
  SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{#eta}", "#it{#varphi}", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.27,512,505);//#it{p}_{T} (GeV/#it{c})
  histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
  // histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
  histo2DAllPi0FWHM->SetMaximum(20);
  histo2DAllPi0FWHM->DrawCopy();


  h_PhiEta_HCal->SetMaximum(20);
  h_PhiEta_HCal->SetLineColor(kBlack);
  h_PhiEta_HCal->SetFillColor(kBlack);
  h_PhiEta_HCal->Draw("box,same");
  h_PhiEta_iEta_HCal->SetLineColor(kBlue+2);
  h_PhiEta_iEta_HCal->SetFillColor(kBlue+2);
  h_PhiEta_iEta_HCal->Draw("box,same");
  h_PhiEta_iPhi_HCal->SetLineColor(kRed+2);
  h_PhiEta_iPhi_HCal->SetFillColor(kRed+2);
  h_PhiEta_iPhi_HCal->Draw("box,same");

  legendPhiEta    = GetAndSetLegend2(0.10, 0.08, 0.38, 0.08+(3*0.85*textsizeLabelsWidth),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendPhiEta->AddEntry(h_PhiEta_ECal, "all towers","plf");
  legendPhiEta->AddEntry(h_PhiEta_iEta_ECal, "towers iEta = 1, 30, 70","fl");
  legendPhiEta->AddEntry(h_PhiEta_iPhi_ECal, "towers iPhi = 1, 50, 100","fl");
  legendPhiEta->Draw();
  
  drawLatexAdd("#it{#bf{ECCE}} simulation",0.10,0.86,textsizeLabelsWidth,kFALSE,kFALSE,false);
  drawLatexAdd(Form("HCals"),0.10,0.78,textsizeLabelsWidth,kFALSE,kFALSE,false);
  drawLatexAdd(Form("a)"),0.97,0.08,textsizeLabelsWidth,kFALSE,kFALSE,true);

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

  TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",800, -4, 4, 400, -TMath::Pi(), TMath::Pi());//, 100.1, 160.9);//125.1, 155.9);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{#eta}", "#it{#varphi}", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.27,512,505);
  histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0Mass->GetXaxis()->SetNoExponent();
  histo2DAllPi0Mass->SetMaximum(20);
  histo2DAllPi0Mass->DrawCopy();


  h_PhiEta_ECal->SetMaximum(20);
  h_PhiEta_ECal->SetLineColor(kBlack);
  h_PhiEta_ECal->SetFillColor(kBlack);
  h_PhiEta_ECal->Draw("box,same");
  h_PhiEta_iEta_ECal->SetLineColor(kBlue+2);
  h_PhiEta_iEta_ECal->SetFillColor(kBlue+2);
  h_PhiEta_iEta_ECal->Draw("box,same");
  h_PhiEta_iPhi_ECal->SetLineColor(kRed+2);
  h_PhiEta_iPhi_ECal->SetFillColor(kRed+2);
  h_PhiEta_iPhi_ECal->Draw("box,same");
  drawLatexAdd(Form("ECals"),0.10,0.88,textsizeLabelsMass,kFALSE,kFALSE,false);
  drawLatexAdd(Form("b)"),0.97,0.19,textsizeLabelsMass,kFALSE,kFALSE,true);

  canvasGeoPaper->Update();
  canvasGeoPaper->Print(Form("%s/GeometryECCE_Acceptance_Paper.%s",outputDir.Data(),suffix.Data()));
}

