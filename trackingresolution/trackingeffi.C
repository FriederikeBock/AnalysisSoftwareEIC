#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void NormalizeByBinWidth(TH1D * hist) {
//	hist->Scale(1,"width");
//	return;
  int nBins = hist->GetNbinsX();
  for (int i = 1; i <= nBins; i++) {
    double fBinWidth = hist->GetXaxis()->GetBinWidth(i);
    hist->SetBinContent(i,hist->GetBinContent(i)/fBinWidth);
    hist->SetBinError(i,hist->GetBinError(i)/fBinWidth);
  }
}

void trackingeffi(
                            TString inputFileName   = "file.root",
                            TString suffix          = "pdf",
                            TString addLabel        = "",
                            Int_t trackCuts         = 0
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  SetPlotStyle();
  TString dateForOutput             = ReturnDateStringForOutput();
  TString readTrackClass            = "";
  TString writeLabel                = "";
  TString labelPlotCuts             = "";
  if (trackCuts == 1){
    readTrackClass                  = "LI2";
    addLabel                        = addLabel+"-LI2";
    writeLabel                      = "LI2";
    labelPlotCuts                   = "#geq 2 tracker hits";
  } else if (trackCuts == 2){
    readTrackClass                  = "LI3";
    addLabel                        = addLabel+"-LI3";
    writeLabel                      = "LI3";
    labelPlotCuts                   = "#geq 3 tracker hits";
  }
    
  for (Int_t eR = 0; eR < 3; eR++) cout << "before: "<< labelEtaRange[eR].Data() << endl;
  
  TString collisionSystem = GetCollisionEnergy(addLabel);
  TString magnetLabel     = GetMagnetLabel(addLabel);
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addLabel.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }

  TString outputDir                 = Form("plots/%s/Effi%s",dateForOutput.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TString detLabel        = GetTrackerLabel(addLabel);
  Int_t nActiveEta        = 14;
  Double_t maxPtSigma     = 0.175;
  Double_t maxEtaSigma    = 0.005;
  Double_t maxPhiSigma    = 0.01;
  if (addLabel.Contains("Hist")){
    maxPtSigma         = 0.42;
    maxEtaSigma        = 0.01;
    maxPhiSigma        = 0.1;
  }
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  TH1D* h_spectra_MC_MCpT[nPID][nEta+1]         = {{NULL}};
  TH1D* h_spectra_rec_pT[nPID][nEta+1]          = {{NULL}};
  TH1D* h_spectra_rec_MCpT[nPID][nEta+1]        = {{NULL}};
  TH1D* h_spectra_MC_MCp[nPID][nEta+1]          = {{NULL}};
  TH1D* h_spectra_rec_p[nPID][nEta+1]           = {{NULL}};
  TH1D* h_spectra_rec_MCp[nPID][nEta+1]         = {{NULL}};
  
  TH1D* h_spectraReb_MC_MCpT[nPID][nEta+1]      = {{NULL}};
  TH1D* h_spectraReb_rec_pT[nPID][nEta+1]       = {{NULL}};
  TH1D* h_spectraReb_rec_MCpT[nPID][nEta+1]     = {{NULL}};
  TH1D* h_spectraReb_MC_MCp[nPID][nEta+1]       = {{NULL}};
  TH1D* h_spectraReb_rec_p[nPID][nEta+1]        = {{NULL}};
  TH1D* h_spectraReb_rec_MCp[nPID][nEta+1]      = {{NULL}};
  
  TH1D* h_effi_rec_pT[nPID][nEta+1]             = {{NULL}};
  TH1D* h_effi_rec_MCpT[nPID][nEta+1]           = {{NULL}};
  TH1D* h_effi_rec_p[nPID][nEta+1]              = {{NULL}};
  TH1D* h_effi_rec_MCp[nPID][nEta+1]            = {{NULL}};
  
  TH2F* h_trackMapMC_eta_MCpT[nPID]             = {NULL};
  TH2F* h_trackMapRec_eta_pT[nPID]              = {NULL};
  TH2F* h_trackMapRec_eta_MCpT[nPID]            = {NULL};
  TH2F* h_trackMapMC_eta_MCp[nPID]              = {NULL};
  TH2F* h_trackMapRec_eta_p[nPID]               = {NULL};
  TH2F* h_trackMapRec_eta_MCp[nPID]             = {NULL};
  
  TFile* inputFile  = new TFile(inputFileName.Data());

  for (Int_t pid = 0; pid < nPID; pid++){
    cout << Form("h_%s_MC_pT", partNameET[pid].Data()) << endl;
    h_trackMapMC_eta_MCpT[pid]                  = (TH2F*)inputFile->Get(Form("h_%s_MC_pT", partNameET[pid].Data()));
    h_trackMapMC_eta_MCpT[pid]->Sumw2();
    cout << Form("h_%s_rec_pT", partNameET[pid].Data()) << endl;
    h_trackMapRec_eta_pT[pid]                   = (TH2F*)inputFile->Get(Form("h_%s_rec_pT", partNameET[pid].Data()));
    h_trackMapRec_eta_pT[pid]->Sumw2();
    cout << Form("h_%s_rec_truepT", partNameET[pid].Data()) << endl;
    h_trackMapRec_eta_MCpT[pid]                 = (TH2F*)inputFile->Get(Form("h_%s_rec_truepT", partNameET[pid].Data()));
    h_trackMapRec_eta_MCpT[pid]->Sumw2();
    cout << Form("h_%s_MC_p", partNameET[pid].Data()) << endl;
    h_trackMapMC_eta_MCp[pid]                   = (TH2F*)inputFile->Get(Form("h_%s_MC_p", partNameET[pid].Data()));
    h_trackMapMC_eta_MCp[pid]->Sumw2();
    cout << Form("h_%s_rec_p", partNameET[pid].Data()) << endl;
    h_trackMapRec_eta_p[pid]                    = (TH2F*)inputFile->Get(Form("h_%s_rec_p", partNameET[pid].Data()));
    h_trackMapRec_eta_p[pid]->Sumw2();
    cout << Form("h_%s_rec_truep", partNameET[pid].Data()) << endl;
    h_trackMapRec_eta_MCp[pid]                  = (TH2F*)inputFile->Get(Form("h_%s_rec_truep", partNameET[pid].Data()));
    h_trackMapRec_eta_MCp[pid]->Sumw2();
  }
  
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput)std::cout << Form("%1.1f<#eta<%1.1f",etaMin,etaMax)  << std::endl;
    
    for (Int_t pid = 0; pid < nPID; pid++){
      h_spectra_MC_MCpT[pid][iEta]         = (TH1D*)h_trackMapMC_eta_MCpT[pid]->ProjectionX(Form("spectraMC%s_MCpT_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapMC_eta_MCpT[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapMC_eta_MCpT[pid]->GetYaxis()->FindBin(etaMax-0.001),"e"); 
      h_spectraReb_MC_MCpT[pid][iEta]      = (TH1D*)h_spectra_MC_MCpT[pid][iEta]->Rebin(nPt, Form("spectraRebMC%s_MCpT_%d",partName[pid].Data(), iEta),
                                                                                          partPt);
      NormalizeByBinWidth(h_spectra_MC_MCpT[pid][iEta]);
      NormalizeByBinWidth(h_spectraReb_MC_MCpT[pid][iEta]);
      DrawGammaSetMarker(h_spectra_MC_MCpT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraReb_MC_MCpT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      
      h_spectra_rec_pT[pid][iEta]          = (TH1D*)h_trackMapRec_eta_pT[pid]->ProjectionX(Form("spectraRec%s_pT_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapRec_eta_pT[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapRec_eta_pT[pid]->GetYaxis()->FindBin(etaMax-0.001),"e"); 
      h_spectraReb_rec_pT[pid][iEta]       = (TH1D*)h_spectra_rec_pT[pid][iEta]->Rebin(nPt, Form("spectraRebRec%s_pT_%d",partName[pid].Data(), iEta),
                                                                                          partPt);
      NormalizeByBinWidth(h_spectra_rec_pT[pid][iEta]);
      NormalizeByBinWidth(h_spectraReb_rec_pT[pid][iEta]);
      DrawGammaSetMarker(h_spectra_rec_pT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraReb_rec_pT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

      h_spectra_rec_MCpT[pid][iEta]        = (TH1D*)h_trackMapRec_eta_MCpT[pid]->ProjectionX(Form("spectraRec%s_MCpT_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapRec_eta_MCpT[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapRec_eta_MCpT[pid]->GetYaxis()->FindBin(etaMax-0.001),"e"); 
      h_spectraReb_rec_MCpT[pid][iEta]     = (TH1D*)h_spectra_rec_MCpT[pid][iEta]->Rebin(nPt, Form("spectraRebRec%s_MCpT_%d",partName[pid].Data(), iEta),
                                                                                          partPt);
      NormalizeByBinWidth(h_spectra_rec_MCpT[pid][iEta]);
      NormalizeByBinWidth(h_spectraReb_rec_MCpT[pid][iEta]);
      DrawGammaSetMarker(h_spectra_rec_MCpT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraReb_rec_MCpT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_spectra_MC_MCp[pid][iEta]          = (TH1D*)h_trackMapMC_eta_MCp[pid]->ProjectionX(Form("spectraMC%s_MCp_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapMC_eta_MCp[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapMC_eta_MCp[pid]->GetYaxis()->FindBin(etaMax-0.001),"e");       
      h_spectraReb_MC_MCp[pid][iEta]       = (TH1D*)h_spectra_MC_MCp[pid][iEta]->Rebin(nP-2, Form("spectraRebMC%s_MCp_%d",partName[pid].Data(), iEta),
                                                                                          partP);
      NormalizeByBinWidth(h_spectra_MC_MCp[pid][iEta]);
      NormalizeByBinWidth(h_spectraReb_MC_MCp[pid][iEta]);
      DrawGammaSetMarker(h_spectra_MC_MCp[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraReb_MC_MCp[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

      h_spectra_rec_p[pid][iEta]           = (TH1D*)h_trackMapRec_eta_p[pid]->ProjectionX(Form("spectraRec%s_p_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapRec_eta_p[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapRec_eta_p[pid]->GetYaxis()->FindBin(etaMax-0.001),"e"); 
      h_spectraReb_rec_p[pid][iEta]        = (TH1D*)h_spectra_rec_p[pid][iEta]->Rebin(nP-2, Form("spectraRebRec%s_p_%d",partName[pid].Data(), iEta),
                                                                                          partP);
      NormalizeByBinWidth(h_spectra_rec_p[pid][iEta]);
      NormalizeByBinWidth(h_spectraReb_rec_p[pid][iEta]);
      DrawGammaSetMarker(h_spectra_rec_p[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraReb_rec_p[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

      h_spectra_rec_MCp[pid][iEta]         = (TH1D*)h_trackMapRec_eta_MCp[pid]->ProjectionX(Form("spectraRec%s_MCp_%d",partName[pid].Data(), iEta), 
                                                                          h_trackMapRec_eta_MCp[pid]->GetYaxis()->FindBin(etaMin+0.001), h_trackMapRec_eta_MCp[pid]->GetYaxis()->FindBin(etaMax-0.001),"e"); 
      h_spectraReb_rec_MCp[pid][iEta]      = (TH1D*)h_spectra_rec_MCp[pid][iEta]->Rebin(nP-2, Form("spectraRebRec%s_MCp_%d",partName[pid].Data(), iEta),
                                                                                          partP);
      NormalizeByBinWidth(h_spectra_rec_MCp[pid][iEta]);
      NormalizeByBinWidth(h_spectraReb_rec_MCp[pid][iEta]);
      DrawGammaSetMarker(h_spectra_rec_MCp[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      DrawGammaSetMarker(h_spectraReb_rec_MCp[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);

      h_effi_rec_pT[pid][iEta]             = (TH1D*)h_spectraReb_rec_pT[pid][iEta]->Clone(Form("effi%s_pT_%d",partName[pid].Data(), iEta));
      h_effi_rec_pT[pid][iEta]->Divide(h_spectraReb_rec_pT[pid][iEta],h_spectraReb_MC_MCpT[pid][iEta],1,1,"B");
      h_effi_rec_MCpT[pid][iEta]           = (TH1D*)h_spectraReb_rec_MCpT[pid][iEta]->Clone(Form("effi%s_MCpT_%d",partName[pid].Data(), iEta));
      h_effi_rec_MCpT[pid][iEta]->Divide(h_spectraReb_rec_MCpT[pid][iEta],h_spectraReb_MC_MCpT[pid][iEta],1,1,"B");
      h_effi_rec_p[pid][iEta]              = (TH1D*)h_spectraReb_rec_p[pid][iEta]->Clone(Form("effi%s_p_%d",partName[pid].Data(), iEta));
      h_effi_rec_p[pid][iEta]->Divide(h_spectraReb_rec_p[pid][iEta],h_spectraReb_MC_MCp[pid][iEta],1,1,"B");
      h_effi_rec_MCp[pid][iEta]            = (TH1D*)h_spectraReb_rec_MCp[pid][iEta]->Clone(Form("effi%s_MCp_%d",partName[pid].Data(), iEta));
      h_effi_rec_MCp[pid][iEta]->Divide(h_spectraReb_rec_MCp[pid][iEta],h_spectraReb_MC_MCp[pid][iEta],1,1,"B");
//       cout << h_effi_rec_pT[pid][iEta] << "\t" << h_effi_rec_MCpT[pid][iEta] << "\t" << h_effi_rec_p[pid][iEta] << "\t" << h_effi_rec_MCp[pid][iEta] << endl;
    }
  }
  
  // 1D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
  DrawGammaCanvasSettings( cReso, 0.085, 0.025, 0.02, 0.105);
  // cReso->SetLogz();

  TH2F* histoDummyEffiPt   = new TH2F("histoDummyEffiPt","histoDummyEffiPt",1000,0, 20,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiPt, "#it{p}_{T} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);
  histoDummyEffiPt->GetXaxis()->SetNoExponent();
  histoDummyEffiPt->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiPt->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyEffiMCPt   = new TH2F("histoDummyEffiMCPt","histoDummyEffiMCPt",1000,0, 20,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiMCPt, "#it{p}_{T}^{MC} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);;
  histoDummyEffiMCPt->GetXaxis()->SetNoExponent();
  histoDummyEffiMCPt->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiMCPt->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyEffiP   = new TH2F("histoDummyEffiP","histoDummyEffiP",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiP, "#it{p} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);
  histoDummyEffiP->GetXaxis()->SetNoExponent();
  histoDummyEffiP->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiP->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyEffiMCP   = new TH2F("histoDummyEffiMCP","histoDummyEffiMCP",1000,0, 100,1000,0.0, 1.35);
  SetStyleHistoTH2ForGraphs(histoDummyEffiMCP, "#it{p}^{MC} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);
  histoDummyEffiMCP->GetXaxis()->SetNoExponent();
  histoDummyEffiMCP->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEffiMCP->GetXaxis()->SetMoreLogLabels(kTRUE);

  TLegend* legendPtResM     = GetAndSetLegend2(0.45, 0.14, 0.95, 0.14+(nActiveEta/3*0.75*textSizeLabelsRel),0.75*textSizeLabelsPixel, 3, "", 43, 0.2);
  TLegend* legendPtResMP    = GetAndSetLegend2(0.45, 0.14, 0.95, 0.14+(nActiveEta/3*0.75*textSizeLabelsRel),0.75*textSizeLabelsPixel, 3, "", 43, 0.2);
  TLegend* legendPtResPID   = GetAndSetLegend2(0.75, 0.14, 0.95, 0.14+(3*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  
  for (Int_t pid = 0; pid < nPID; pid++){
    cout << "PID: "<< pid << endl;
    histoDummyEffiPt->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_effi_rec_pT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_effi_rec_pT[pid][iEta]->Draw("same,p");
      if (pid == 0){
        if (iEta == nEta )
          legendPtResM->AddEntry(h_effi_rec_pT[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
        else 
          legendPtResM->AddEntry(h_effi_rec_pT[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
    }
    legendPtResM->Draw();
    DrawGammaLines(0, 20, 1., 1., 2, kGray+2, 7);
      
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/Effi_%s_pT.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));

    histoDummyEffiMCPt->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_effi_rec_MCpT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_effi_rec_MCpT[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    DrawGammaLines(0, 20, 1., 1., 2, kGray+2, 7);
      
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/Effi_%s_MCpT.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));

    cReso->SetLogx();
    histoDummyEffiP->Draw();
    DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_effi_rec_p[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_effi_rec_p[pid][iEta]->Draw("same,p");
      if (pid == 0){
        if (iEta == nEta )
          legendPtResMP->AddEntry(h_effi_rec_p[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
        else 
          legendPtResMP->AddEntry(h_effi_rec_p[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
    }
    legendPtResMP->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/Effi_%s_p.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));

    histoDummyEffiMCP->Draw();
    DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_effi_rec_MCp[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_effi_rec_MCp[pid][iEta]->Draw("same,p");
    }
    legendPtResMP->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    cReso->Print(Form("%s/Effi_%s_MCp.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    cReso->SetLogx(kFALSE);
  }
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    histoDummyEffiMCPt->Draw();
    legendPtResPID->Clear();
    DrawGammaLines(0, 20, 1., 1., 2, kGray+2, 7);
    for (Int_t pid =0; pid < nPID; pid++){
      
      DrawGammaSetMarker(h_effi_rec_MCpT[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_effi_rec_MCpT[pid][iEta]->Draw("same,p");      
      legendPtResPID->AddEntry(h_effi_rec_MCpT[pid][iEta],partLabel[pid].Data(),"p");
    }
    legendPtResPID->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s", etaMin, etaMax, detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/EffiPID_MCpT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    cReso->SetLogx();
    histoDummyEffiMCP->Draw();
    DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
    for (Int_t pid =0; pid < nPID; pid++){
      
      DrawGammaSetMarker(h_effi_rec_MCp[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_effi_rec_MCp[pid][iEta]->Draw("same,p");      
    }
    legendPtResPID->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s", etaMin, etaMax, detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/EffiPID_MCp_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    cReso->SetLogx(kFALSE);
  }

  
  TLegend* legendPEtaRegion[3]  = {NULL, NULL, NULL}; 
  TLegend* legendPtEtaRegion[3]  = {NULL, NULL, NULL}; 
  cReso->SetLogx();
  for (Int_t pid = 0; pid < nPID; pid++){
    
    for (Int_t eR = 0; eR < 3; eR++){
      histoDummyEffiMCPt->Draw();
//       cout << labelEtaRange[eR].Data() << endl;
      if (pid == 0) legendPtEtaRegion[eR]  = GetAndSetLegend2(0.75, 0.14, 0.95, 0.14+((maxNEtaBins[eR]+1)*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, labelEtaRange[eR], 43, 0.2);
      for(Int_t iEta=minEtaBin[eR]; iEta<maxEtaBin[eR]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_rec_MCpT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_rec_MCpT[pid][iEta]->Draw("same,p");
        if (pid == 0){
          if (iEta == nEta )
            legendPtEtaRegion[eR]->AddEntry(h_effi_rec_MCpT[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
          else 
            legendPtEtaRegion[eR]->AddEntry(h_effi_rec_MCpT[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
        }
      }
      legendPtEtaRegion[eR]->Draw();
      DrawGammaLines(0., 20, 1., 1., 2, kGray+2, 7);
        
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/Effi_%s_%s_MCpT.%s", outputDir.Data(), nameOutEtaRange[eR].Data(), partName[pid].Data(), suffix.Data()));

      histoDummyEffiPt->Draw();
//       cout << labelEtaRange[eR].Data() << endl;
      for(Int_t iEta=minEtaBin[eR]; iEta<maxEtaBin[eR]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_rec_pT[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_rec_pT[pid][iEta]->Draw("same,p");
      }
      legendPtEtaRegion[eR]->Draw();
      DrawGammaLines(0., 20, 1., 1., 2, kGray+2, 7);
        
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/Effi_%s_%s_pT.%s", outputDir.Data(), nameOutEtaRange[eR].Data(), partName[pid].Data(), suffix.Data()));

      cReso->SetLogx();
      histoDummyEffiMCP->Draw();
//       cout << labelEtaRange[eR].Data() << endl;
      if (pid == 0) legendPEtaRegion[eR]  = GetAndSetLegend2(0.75, 0.14, 0.95, 0.14+((maxNEtaBins[eR]+1)*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, labelEtaRange[eR], 43, 0.2);
      for(Int_t iEta=minEtaBin[eR]; iEta<maxEtaBin[eR]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_rec_MCp[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_rec_MCp[pid][iEta]->Draw("same,p");
        if (pid == 0){
          if (iEta == nEta )
            legendPEtaRegion[eR]->AddEntry(h_effi_rec_MCp[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
          else 
            legendPEtaRegion[eR]->AddEntry(h_effi_rec_MCp[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
        }
      }
      legendPEtaRegion[eR]->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/Effi_%s_%s_MCp.%s", outputDir.Data(), nameOutEtaRange[eR].Data(), partName[pid].Data(), suffix.Data()));

      histoDummyEffiP->Draw();
//       cout << labelEtaRange[eR].Data() << endl;
      for(Int_t iEta=minEtaBin[eR]; iEta<maxEtaBin[eR]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_effi_rec_p[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_effi_rec_p[pid][iEta]->Draw("same,p");
      }
      legendPEtaRegion[eR]->Draw();
      DrawGammaLines(0.1, 100, 1., 1., 2, kGray+2, 7);
        
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/Effi_%s_%s_p.%s", outputDir.Data(), nameOutEtaRange[eR].Data(), partName[pid].Data(), suffix.Data()));

      cReso->SetLogx(kFALSE);
    }
  }
   
  TFile* outputFile  = new TFile(inputFileName.Data(),"UPDATE");
  TDirectoryFile* directoryTrEffi = (TDirectoryFile*)outputFile->Get("TrackingEfficiencyProjections");
  if (!directoryTrEffi){
    outputFile->mkdir("TrackingEfficiencyProjections");
    directoryTrEffi = (TDirectoryFile*)outputFile->Get("TrackingEfficiencyProjections");
  }
  directoryTrEffi->cd();
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    for (Int_t pid = 0; pid < 6; pid++){
      h_effi_rec_pT[pid][iEta]->Write(Form("effi%s_pT_%d",partName[pid].Data(), iEta),TObject::kOverwrite);
      h_effi_rec_MCpT[pid][iEta]->Write(Form("effi%s_MCpT_%d",partName[pid].Data(), iEta),TObject::kOverwrite);
      h_effi_rec_p[pid][iEta]->Write(Form("effi%s_p_%d",partName[pid].Data(), iEta),TObject::kOverwrite);
      h_effi_rec_MCp[pid][iEta]->Write(Form("effi%s_MCp_%d",partName[pid].Data(), iEta),TObject::kOverwrite);
      h_spectra_MC_MCpT[pid][iEta]->Write(Form("spectraMC%s_MCpT_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectra_rec_pT[pid][iEta]->Write(Form("spectraRec%s_pT_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectra_rec_MCpT[pid][iEta]->Write(Form("spectraRec%s_MCpT_%d",partName[pid].Data(), iEta), TObject::kOverwrite); 
      h_spectra_MC_MCp[pid][iEta]->Write( Form("spectraMC%s_MCp_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectra_rec_p[pid][iEta]->Write( Form("spectraRec%s_p_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectra_rec_MCp[pid][iEta]->Write( Form("spectraRec%s_MCp_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectraReb_MC_MCpT[pid][iEta]->Write( Form("spectraRebMC%s_MCpT_%d", partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectraReb_rec_pT[pid][iEta]->Write( Form("spectraRebRec%s_pT_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectraReb_rec_MCpT[pid][iEta]->Write( Form("spectraRebRec%s_MCpT_%d",partName[pid].Data(), iEta), TObject::kOverwrite);      
      h_spectraReb_MC_MCp[pid][iEta]->Write( Form("spectraRebMC%s_MCp_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectraReb_rec_p[pid][iEta]->Write( Form("spectraRebRec%s_p_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
      h_spectraReb_rec_MCp[pid][iEta]->Write( Form("spectraRebRec%s_MCp_%d",partName[pid].Data(), iEta), TObject::kOverwrite);
    }
  }
  outputFile->Write();
  outputFile->Close();
}
