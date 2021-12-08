#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 

//************************************************************************************
//****************** Main function ***************************************************
//************************************************************************************
void showershapeStudies(
    TString inputFileNameE    = "central_allculsters.root",
    TString inputFileNamePi   = "central_allculsters.root",
    TString suffix            = "pdf",
    Bool_t doFitting          = kTRUE,
    TString collisionsSys     = "PythiaMB"
){  
  TString clusterizerName = "MA clusters";
  const int maxcalo = 3;
  TString caloName[maxcalo+1] = {"EEMC", "BECAL", "FEMC"};
  TString caloNamePlot[maxcalo+1] = {"EEMC", "BEMC", "FEMC"};
  bool caloactive[maxcalo] = {false};
  Double_t mass_pi0PDG = 0.1349770;
  double eopcutvalue[maxcalo] = {0.8, 0.75, 0.77};
  double eopcutvaluemax[maxcalo] = {1.9, 1.9, 1.9};
  double etarangemin[maxcalo] = {-4, -1.7, 1.3};
  double etarangemax[maxcalo] = {-1.7, 1.3, 4.0};
  double energymax[maxcalo] = {50, 20, 50.0};
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem           = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString perfLabel                 = "#it{#bf{ECCE}} simulation";
  TString dateForOutput             = ReturnDateStringForOutput();

  Color_t colorCalo[maxcalo];
  colorCalo[0] = getCaloColor("EEMC", 0);
  colorCalo[1] = getCaloColor("BEMC", 0);
  colorCalo[2] = getCaloColor("FEMC", 0);
  
  Color_t colorCalo_light[maxcalo];
  colorCalo_light[0] = getCaloColor("EEMC", 1);
  colorCalo_light[1] = getCaloColor("BEMC", 1);
  colorCalo_light[2] = getCaloColor("FEMC", 1);
  
  Color_t colorCalo_light2[nPID]         = {kBlue-2, kRed-2, kGreen-2, kCyan-6, kBlue+1, kOrange};

  if (collisionsSys.CompareTo("PythiaMB") == 0){
    collisionSystem   = "e-p: 18#times 275 GeV";
  } else if (collisionsSys.CompareTo("SingleElectron") == 0){
    collisionSystem   = "e^{-}";
  } else if (collisionsSys.CompareTo("SingleKaon") == 0){
    collisionSystem   = "K^{-}";
  } else if (collisionsSys.CompareTo("SinglePion") == 0){
    collisionSystem   = "#pi^{-}";
  } else if (collisionsSys.CompareTo("SingleProton") == 0){
    collisionSystem   = "p";
  } else if (collisionsSys.CompareTo("SingleNeutron") == 0){
    collisionSystem   = "n";
  } else if (collisionsSys.CompareTo("SinglePart") == 0){
    collisionSystem   = "single particles";
  }

  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;

  const int Nusefilltrue = 2;
  TString str_EoverP_FillTrue[Nusefilltrue] = {"recP","trueP"};
  
  TString outputDir                 = Form("plotsshowershape/%s",dateForOutput.Data());
  // if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFile                      = new TFile(inputFileNameE.Data());
  if(inputFile->IsZombie()){
    cout << "inputFile not found!" << endl;
  }
  TFile* inputFilePi                      = new TFile(inputFileNamePi.Data());
  if(inputFilePi->IsZombie()){
    cout << "inputFilePi not found!" << endl;
  }

  //************************** Read data **************************************************
  const int nEne = 16;
  const static Double_t partE[]   = {   0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2, 3, 4, 5, 6, 8, 10,15,20,25,  30};

  int padebin = 0;
  int padebinelec = 0;
  const int maxmatch = 2;
  TH2F*  h_M02_E_pions[maxcalo]      = {NULL}; 
  TH2F*  h_M02_E_electrons[maxcalo]  = {NULL}; 
  TH2F*  h_M20_E_pions[maxcalo]      = {NULL}; 
  TH2F*  h_M20_E_electrons[maxcalo]  = {NULL}; 

  TH1F* histM02_proj_pions[maxcalo][50]     = {{NULL}};
  TH1F* histM02_proj_electrons[maxcalo][50] = {{NULL}};
  TH1F* histM20_proj_pions[maxcalo][50]     = {{NULL}};
  TH1F* histM20_proj_electrons[maxcalo][50] = {{NULL}};

  double M02_cut_90[maxcalo][50]        = {{0}};
  double pirejectM02_90[maxcalo][50]    = {{0}};
  double M20_cut_90[maxcalo][50]        = {{0}};
  double pirejectM20_90[maxcalo][50]    = {{0}};

  
  TH2F * histProjDummy                           = new TH2F("histProjDummy","histProjDummy",1000,0, 4,1000,0.8, 1e6);
  SetStyleHistoTH2ForGraphs(histProjDummy, "#sigma_{long}^{2}", "counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  TH2F * histProjDummySA                           = new TH2F("histProjDummySA","histProjDummySA",1000,0, 4,1000,0.8, 1e6);
  SetStyleHistoTH2ForGraphs(histProjDummySA, "#sigma_{short}^{2}", "counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  // TString matchbasis = "trackbased";
  int imatch=1;
  for(int icalo=0;icalo<maxcalo;icalo++){

    gSystem->Exec(Form("mkdir -p %s/%s",outputDir.Data(), caloNamePlot[icalo].Data()));
    gSystem->Exec(Form("mkdir -p %s/%s",outputDir.Data(), caloNamePlot[icalo].Data()));
    cout << "processing " << caloName[icalo].Data() << endl;
    
    TString NameM02E = Form("%s/h_CS_clusters_M02_part_E_truth_%s_MA_electron",caloName[icalo].Data(), caloName[icalo].Data() );
    h_M02_E_electrons[icalo]    = (TH2F*)inputFile->Get(NameM02E.Data());
    if(!h_M02_E_electrons[icalo]){
      cout << "\t" <<NameM02E.Data() << " not found" << endl;
      continue;
    } else {
      caloactive[icalo]=true;
    };
    TString NameM02Pi = Form("%s/h_CS_clusters_M02_part_E_truth_%s_MA_cpion",caloName[icalo].Data(), caloName[icalo].Data() );
    h_M02_E_pions[icalo]        = (TH2F*)inputFilePi->Get(NameM02Pi.Data());
    if(!h_M02_E_pions[icalo]){
      cout << "\t" <<NameM02Pi.Data() << " not found" << endl;
      caloactive[icalo]=false;
      continue;
    } else {
      caloactive[icalo]=true;
    };

    TString NameM20E = Form("%s/h_CS_clusters_M20_part_E_truth_%s_MA_electron",caloName[icalo].Data(), caloName[icalo].Data() );
    h_M20_E_electrons[icalo]    = (TH2F*)inputFile->Get(NameM20E.Data());
    TString NameM20Pi = Form("%s/h_CS_clusters_M20_part_E_truth_%s_MA_cpion",caloName[icalo].Data(), caloName[icalo].Data() );
    h_M20_E_pions[icalo]        = (TH2F*)inputFilePi->Get(NameM20Pi.Data());
    
    
    float latextopedge = 0.90;
    float latexsideedge = 0.92;
    float  textSizeLabelsRel      = 55./1200;
    if (h_M02_E_electrons[icalo]->GetEntries() > 100 || h_M02_E_pions[icalo]->GetEntries()>100){
      TCanvas* canvasEoverP2D = new TCanvas(Form("canvasEoverP2D%s",caloName[icalo].Data()),"",0,0,1500,800);
      DrawGammaCanvasSettings( canvasEoverP2D, 0.095, 0.01, 0.01, 0.105);
      canvasEoverP2D->Divide(2,1);
      canvasEoverP2D->cd(1);
      canvasEoverP2D->cd(1)->SetLogy();
      canvasEoverP2D->cd(1)->SetLogz();
      h_M02_E_electrons[icalo]->GetYaxis()->SetRangeUser(0.1,30);
      SetStyleHistoTH2ForGraphs(h_M02_E_electrons[icalo],"#sigma_{long}^{2}", "E_{MC} GeV", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      h_M02_E_electrons[icalo]->Draw("colz");
      drawLatexAdd("true e^{-}",0.17,latextopedge,textSizeLabelsRel,false,false,false);      
      canvasEoverP2D->cd(2);
      canvasEoverP2D->cd(2)->SetLogy();
      canvasEoverP2D->cd(2)->SetLogz();
      
      h_M02_E_pions[icalo]->GetYaxis()->SetRangeUser(0.1,30);
      SetStyleHistoTH2ForGraphs(h_M02_E_pions[icalo],"#sigma_{long}^{2}", "E_{MC} GeV", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      h_M02_E_pions[icalo]->Draw("colz");
    
      drawLatexAdd(perfLabel.Data(),latexsideedge,latextopedge,textSizeLabelsRel,false,false,true);
      drawLatexAdd(collisionSystem.Data(),latexsideedge,latextopedge-0.05,textSizeLabelsRel,false,false,true);
      drawLatexAdd(Form("Calorimeter: %s", caloNamePlot[icalo].Data()),latexsideedge,latextopedge-0.1,textSizeLabelsRel,false,false,true);      
      drawLatexAdd("true #pi^{-}",0.17,latextopedge,textSizeLabelsRel,false,false,false);      
      canvasEoverP2D->Print(Form("%s/%s/M02_vs_E_%s_%s_trueP.%s", outputDir.Data(), caloNamePlot[icalo].Data(), caloNamePlot[icalo].Data(), clusterizerName.Data(), suffix.Data()));
    }
  if (h_M20_E_electrons[icalo]->GetEntries() > 100 || h_M20_E_pions[icalo]->GetEntries()>100){
      TCanvas* canvasEoverP2D = new TCanvas(Form("canvasEoverP2D%s",caloName[icalo].Data()),"",0,0,1500,800);
      DrawGammaCanvasSettings( canvasEoverP2D, 0.095, 0.01, 0.01, 0.105);
      canvasEoverP2D->Divide(2,1);
      canvasEoverP2D->cd(1);
      canvasEoverP2D->cd(1)->SetLogy();
      canvasEoverP2D->cd(1)->SetLogz();
      h_M20_E_electrons[icalo]->GetYaxis()->SetRangeUser(0.1,30);
      SetStyleHistoTH2ForGraphs(h_M20_E_electrons[icalo],"#sigma_{short}^{2}", "E_{MC} GeV", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      h_M20_E_electrons[icalo]->Draw("colz");
      drawLatexAdd("true e^{-}",0.17,latextopedge,textSizeLabelsRel,false,false,false);      
      canvasEoverP2D->cd(2);
      canvasEoverP2D->cd(2)->SetLogy();
      canvasEoverP2D->cd(2)->SetLogz();
      
      h_M20_E_pions[icalo]->GetYaxis()->SetRangeUser(0.1,30);
      SetStyleHistoTH2ForGraphs(h_M20_E_pions[icalo],"#sigma_{short}^{2}", "E_{MC} GeV", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      h_M20_E_pions[icalo]->Draw("colz");
    
      drawLatexAdd(perfLabel.Data(),latexsideedge,latextopedge,textSizeLabelsRel,false,false,true);
      drawLatexAdd(collisionSystem.Data(),latexsideedge,latextopedge-0.05,textSizeLabelsRel,false,false,true);
      drawLatexAdd(Form("Calorimeter: %s", caloNamePlot[icalo].Data()),latexsideedge,latextopedge-0.1,textSizeLabelsRel,false,false,true);      
      drawLatexAdd("true #pi^{-}",0.17,latextopedge,textSizeLabelsRel,false,false,false);      
      canvasEoverP2D->Print(Form("%s/%s/M20_vs_E_%s_%s_trueP.%s", outputDir.Data(), caloNamePlot[icalo].Data(), caloNamePlot[icalo].Data(), clusterizerName.Data(), suffix.Data()));
    }

    
    for(int ebin=0;ebin<nEne-1;ebin++){
      histM02_proj_pions[icalo][ebin]    = (TH1F*)h_M02_E_pions[icalo]->ProjectionX(Form("Pi_%1.1f_%s_x",partE[ebin],caloName[icalo].Data()),h_M02_E_pions[icalo]->GetYaxis()->FindBin(partE[ebin]-0.001), h_M02_E_electrons[icalo]->GetYaxis()->FindBin(partE[ebin+1]+0.001),"e");
      histM02_proj_electrons[icalo][ebin]    = (TH1F*)h_M02_E_electrons[icalo]->ProjectionX(Form("E_%1.1f_%s_x",partE[ebin],caloName[icalo].Data()),h_M02_E_pions[icalo]->GetYaxis()->FindBin(partE[ebin]-0.001), h_M02_E_electrons[icalo]->GetYaxis()->FindBin(partE[ebin+1]+0.001),"e");
      
      Double_t intelectronsall  = histM02_proj_electrons[icalo][ebin]->GetEntries();
      Double_t intpionsall      = histM02_proj_pions[icalo][ebin]->GetEntries();
      Double_t intelectrons_90 = 0;
      M02_cut_90[icalo][ebin] = 0;
      for(int iwindow=1;iwindow<histM02_proj_electrons[icalo][ebin]->GetNbinsX();iwindow++){
        if(abs(intelectrons_90/intelectronsall)>0.90) break;
        else {
          intelectrons_90 = histM02_proj_electrons[icalo][ebin]->Integral(1,iwindow);
          M02_cut_90[icalo][ebin] = histM02_proj_electrons[icalo][ebin]->GetBinCenter(iwindow)+0.5*histM02_proj_electrons[icalo][ebin]->GetBinWidth(iwindow);
        }
      }
      Double_t intpions90e      = histM02_proj_pions[icalo][ebin]->Integral(histM02_proj_pions[icalo][ebin]->GetXaxis()->FindBin(M02_cut_90[icalo][ebin]),histM02_proj_pions[icalo][ebin]->GetNbinsX());
      pirejectM02_90[icalo][ebin] = intpions90e/intpionsall;
      
      histM20_proj_pions[icalo][ebin]    = (TH1F*)h_M20_E_pions[icalo]->ProjectionX(Form("Pi_M20_%1.1f_%s_x",partE[ebin],caloName[icalo].Data()),h_M20_E_pions[icalo]->GetYaxis()->FindBin(partE[ebin]-0.001), h_M20_E_electrons[icalo]->GetYaxis()->FindBin(partE[ebin+1]+0.001),"e");
      histM20_proj_electrons[icalo][ebin]    = (TH1F*)h_M20_E_electrons[icalo]->ProjectionX(Form("E_M20_%1.1f_%s_x",partE[ebin],caloName[icalo].Data()),h_M20_E_pions[icalo]->GetYaxis()->FindBin(partE[ebin]-0.001), h_M20_E_electrons[icalo]->GetYaxis()->FindBin(partE[ebin+1]+0.001),"e");
      
      Double_t intelectronsSA_90 = 0;
      M20_cut_90[icalo][ebin] = 0;
      for(int iwindow=1;iwindow<histM20_proj_electrons[icalo][ebin]->GetNbinsX();iwindow++){
        if(abs(intelectronsSA_90/intelectronsall)>0.90) break;
        else {
          intelectronsSA_90 = histM20_proj_electrons[icalo][ebin]->Integral(1,iwindow);
          M20_cut_90[icalo][ebin] = histM20_proj_electrons[icalo][ebin]->GetBinCenter(iwindow)+0.5*histM20_proj_electrons[icalo][ebin]->GetBinWidth(iwindow);
        }
      }
      Double_t intpions90eSA      = histM20_proj_pions[icalo][ebin]->Integral(histM20_proj_pions[icalo][ebin]->GetXaxis()->FindBin(M20_cut_90[icalo][ebin]),histM20_proj_pions[icalo][ebin]->GetNbinsX());
      pirejectM20_90[icalo][ebin] = intpions90eSA/intpionsall;
    }
    
    float minEtaProj = nominalEtaRegion[icalo][0];
    float maxEtaProj = nominalEtaRegion[icalo][1];

    TCanvas* canvasSSSingleBins;
    split_canvas(canvasSSSingleBins, "canvasSSSingleBins", nEne+1);
    int padnumber = 0;
    for(int z=0; z< nEne ; z++ ){
      if(z==0){
        canvasSSSingleBins->cd(z+1);
        canvasSSSingleBins->cd(z+1)->Clear();
        drawLatexAdd(perfLabel,0.03,0.9,2*textSizeSinglePad,kFALSE);
        drawLatexAdd(Form("Calorimeter: %s ", caloNamePlot[icalo].Data()),0.03,0.9-2*textSizeSinglePad,2*textSizeSinglePad,kFALSE);
        drawLatexAdd(Form("Clusterization type: %s ", clusterizerName.Data()),0.03,0.9-2*2*textSizeSinglePad,2*textSizeSinglePad,kFALSE);
        drawLatexAdd(Form("%1.2f < #eta < %1.1f", minEtaProj, maxEtaProj),0.03,0.9-3*2*textSizeSinglePad,2*textSizeSinglePad,kFALSE);
        
        TLegend* legend1 = GetAndSetLegend2(0.03, 0.1, 0.5, 0.1+2*2*textSizeSinglePad,2*textSizeSinglePad,1,"",42); 
        DrawGammaSetMarker(histM02_proj_electrons[icalo][padnumber],33,1, kBlack, kBlack);
        DrawGammaSetMarker(histM02_proj_pions[icalo][padnumber],28,1, kRed+2, kRed+2);
        legend1->AddEntry(histM02_proj_electrons[icalo][padnumber], "e^{-}", "p");
        legend1->AddEntry(histM02_proj_pions[icalo][padnumber], "#pi^{-}", "p");
        legend1->Draw();
      } else{
        canvasSSSingleBins->cd(z+1);
        DrawVirtualPadSettings( canvasSSSingleBins->cd(z+1), 0.1, 0.02, 0.015, 0.115);
        gPad->SetLogy();
        cout << histM02_proj_electrons[icalo][padnumber]->GetEntries() << endl;
        histProjDummy->GetYaxis()->SetRangeUser(0.8,histM02_proj_electrons[icalo][padnumber]->GetMaximum()*10);
        histProjDummy->DrawCopy();
        if (histM02_proj_electrons[icalo][padnumber]){
          DrawGammaSetMarker(histM02_proj_electrons[icalo][padnumber],33,2, kBlack, kBlack);
          histM02_proj_electrons[icalo][padnumber]->Draw("same,ep");
        }
        if(histM02_proj_pions[icalo][padnumber]){
          DrawGammaSetMarker(histM02_proj_pions[icalo][padnumber],28,1, kRed+2, kRed+2);
          histM02_proj_pions[icalo][padnumber]->Draw("ep,same");
        }

        DrawGammaLines(M02_cut_90[icalo][padnumber], M02_cut_90[icalo][padnumber], 0, 0.5*histM02_proj_electrons[icalo][padnumber]->GetMaximum(),1, kGray+1,3);
        drawLatexAdd(Form("%1.1f < #it{p}_{track} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        padnumber++;
      }
    }
    canvasSSSingleBins->Print(Form("%s/%s/M02_bin_electrons_and_pions_%s_%s_full_trueP.%s", outputDir.Data(), caloNamePlot[icalo].Data(), caloNamePlot[icalo].Data(), clusterizerName.Data(), suffix.Data()));
    
    padnumber = 0;
    for(int z=0; z< nEne ; z++ ){
      if(z==0){
        canvasSSSingleBins->cd(z+1);
        canvasSSSingleBins->cd(z+1)->Clear();
        drawLatexAdd(perfLabel,0.03,0.9,2*textSizeSinglePad,kFALSE);
        drawLatexAdd(Form("Calorimeter: %s ", caloNamePlot[icalo].Data()),0.03,0.9-2*textSizeSinglePad,2*textSizeSinglePad,kFALSE);
        drawLatexAdd(Form("Clusterization type: %s ", clusterizerName.Data()),0.03,0.9-2*2*textSizeSinglePad,2*textSizeSinglePad,kFALSE);
        drawLatexAdd(Form("%1.2f < #eta < %1.1f", minEtaProj, maxEtaProj),0.03,0.9-3*2*textSizeSinglePad,2*textSizeSinglePad,kFALSE);
        
        TLegend* legend1 = GetAndSetLegend2(0.03, 0.1, 0.5, 0.1+2*2*textSizeSinglePad,2*textSizeSinglePad,1,"",42); 
        DrawGammaSetMarker(histM20_proj_electrons[icalo][padnumber],33,1, kBlack, kBlack);
        DrawGammaSetMarker(histM20_proj_pions[icalo][padnumber],28,1, kRed+2, kRed+2);
        legend1->AddEntry(histM20_proj_electrons[icalo][padnumber], "e^{-}", "p");
        legend1->AddEntry(histM20_proj_pions[icalo][padnumber], "#pi^{-}", "p");
        legend1->Draw();
      } else{
        canvasSSSingleBins->cd(z+1);
        DrawVirtualPadSettings( canvasSSSingleBins->cd(z+1), 0.1, 0.02, 0.015, 0.115);
        gPad->SetLogy();
        cout << histM20_proj_electrons[icalo][padnumber]->GetEntries() << endl;
        histProjDummySA->GetYaxis()->SetRangeUser(0.8,histM20_proj_electrons[icalo][padnumber]->GetMaximum()*10);
        histProjDummySA->DrawCopy();
        if (histM20_proj_electrons[icalo][padnumber]){
          DrawGammaSetMarker(histM20_proj_electrons[icalo][padnumber],33,2, kBlack, kBlack);
          histM20_proj_electrons[icalo][padnumber]->Draw("same,ep");
        }
        if(histM20_proj_pions[icalo][padnumber]){
          DrawGammaSetMarker(histM20_proj_pions[icalo][padnumber],28,1, kRed+2, kRed+2);
          histM20_proj_pions[icalo][padnumber]->Draw("ep,same");
        }

        DrawGammaLines(M20_cut_90[icalo][padnumber], M20_cut_90[icalo][padnumber], 0, 0.5*histM20_proj_electrons[icalo][padnumber]->GetMaximum(),1, kGray+1,3);
        drawLatexAdd(Form("%1.1f < #it{p}_{track} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        padnumber++;
      }
    }
    canvasSSSingleBins->Print(Form("%s/%s/M20_bin_electrons_and_pions_%s_%s_full_trueP.%s", outputDir.Data(), caloNamePlot[icalo].Data(), caloNamePlot[icalo].Data(), clusterizerName.Data(), suffix.Data()));
  }
  
  
  
  Style_t markerStyleCalo[maxcalo] = {20, 47, 33};
  Style_t markerStyleCalo2[maxcalo] = {24, 46, 27};
  Style_t markerStyleCalo3[maxcalo] = {25, 42, 28};
  Size_t markerSizeCalo[maxcalo] = {1.5, 1.8, 2.2};
  float toplegval = 0.19;
    
  TCanvas* canvasSScutValue       = new TCanvas("canvasSScutValue", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSScutValue,  0.1, 0.01, 0.015, 0.095);
  canvasSScutValue->SetLogx(1);
    TH2F * histoM02cutValue         = new TH2F("histoM02cutValue", "histoM02cutValue",1000, 0.51, 39, 1000, 0., 2.01 );
    SetStyleHistoTH2ForGraphs( histoM02cutValue, "#it{E}_{true} (GeV/#it{c})", "#it{#sigma}_{long}^{2} cut value",
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
    histoM02cutValue->GetYaxis()->SetLabelOffset(0.001);
    histoM02cutValue->GetXaxis()->SetNoExponent();
    histoM02cutValue->GetXaxis()->SetMoreLogLabels(kTRUE);
    TGraphErrors* graphM02cutValue_90[maxcalo]={NULL};

    histoM02cutValue->DrawCopy();
    TLegend* legendM02Cut3      = GetAndSetLegend2(0.70, 0.95-textSizeLabelsRel, 0.90, 0.95-(textSizeLabelsRel +3*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,1,"",42); 
    TLegend* legendM02Cut4      = GetAndSetLegend2(0.75, 0.95-textSizeLabelsRel, 0.95, 0.95-(textSizeLabelsRel +3*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,1,"",42); 
    for (Int_t icalo = 0; icalo < maxcalo; icalo++){
      graphM02cutValue_90[icalo]   = new TGraphErrors(nEne,partE, M02_cut_90[icalo]);
      while(graphM02cutValue_90[icalo]->GetN()>0 && graphM02cutValue_90[icalo]->GetY()[graphM02cutValue_90[icalo]->GetN()-1] <= 0) graphM02cutValue_90[icalo]->RemovePoint(graphM02cutValue_90[icalo]->GetN()-1);
      DrawGammaSetMarkerTGraph(graphM02cutValue_90[icalo], markerStyleCalo3[icalo], 1.5*markerSizeCalo[icalo], colorCalo_light2[icalo] , colorCalo_light2[icalo]);
      graphM02cutValue_90[icalo]->SetLineStyle(icalo+1);
      graphM02cutValue_90[icalo]->SetLineWidth(3);
      graphM02cutValue_90[icalo]->Draw("samecx0");

      cout << " done l:" << __LINE__ << endl;
      legendM02Cut3->AddEntry(graphM02cutValue_90[icalo]," ","l");
      legendM02Cut4->AddEntry((TObject*)0,Form("%s",caloNamePlot[icalo].Data()),"");
    }

    legendM02Cut3->Draw();
    legendM02Cut4->Draw();

    drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
    drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
    drawLatexAdd("#epsilon_{90%}",0.70,0.95-0.5*textSizeLabelsRel,textSizeLabelsRel,false,false,false);

    DrawGammaLines(0.51, 39, 1,1,2,kGray+1, 2);
  canvasSScutValue->Update();
  canvasSScutValue->Print(Form("%s/M02_cutVlue_AllCalo_true.%s",outputDir.Data(),suffix.Data()));

  toplegval = 0.90;
    TH2F * histoM20cutValue         = new TH2F("histoM20cutValue", "histoM20cutValue",1000, 0.51, 39, 1000, 0., 1.31 );
    SetStyleHistoTH2ForGraphs( histoM20cutValue, "#it{E}_{true} (GeV/#it{c})", "#it{#sigma}_{short}^{2} cut value",
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
    histoM20cutValue->GetYaxis()->SetLabelOffset(0.001);
    histoM20cutValue->GetXaxis()->SetNoExponent();
    histoM20cutValue->GetXaxis()->SetMoreLogLabels(kTRUE);
    TGraphErrors* graphM20cutValue_90[maxcalo]={NULL};
    histoM20cutValue->DrawCopy();
    
    for (Int_t icalo = 0; icalo < maxcalo; icalo++){    
      graphM20cutValue_90[icalo]   = new TGraphErrors(nEne,partE, M20_cut_90[icalo]);
      while(graphM20cutValue_90[icalo]->GetN()>0 && graphM20cutValue_90[icalo]->GetY()[graphM20cutValue_90[icalo]->GetN()-1] <= 0) graphM20cutValue_90[icalo]->RemovePoint(graphM20cutValue_90[icalo]->GetN()-1);
      DrawGammaSetMarkerTGraph(graphM20cutValue_90[icalo], markerStyleCalo3[icalo], 1.5*markerSizeCalo[icalo], colorCalo_light2[icalo] , colorCalo_light2[icalo]);
      graphM20cutValue_90[icalo]->SetLineStyle(icalo+1);
      graphM20cutValue_90[icalo]->SetLineWidth(3);
      graphM20cutValue_90[icalo]->Draw("samecx0");
    }

    legendM02Cut3->Draw();
    legendM02Cut4->Draw();

    drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
    drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
    drawLatexAdd("#epsilon_{90%}",0.70,0.95-0.5*textSizeLabelsRel,textSizeLabelsRel,false,false,false);

    DrawGammaLines(0.51, 39, 1,1,2,kGray+1, 2);
  canvasSScutValue->Update();
  canvasSScutValue->Print(Form("%s/M20_cutVlue_AllCalo_true.%s",outputDir.Data(),suffix.Data()));
  
  toplegval = 0.90;
  TCanvas* canvasReject       = new TCanvas("canvasReject", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasReject,  0.1, 0.01, 0.015, 0.095);
  canvasReject->SetLogx(1);
    TH2F * histoRejectPi         = new TH2F("histoRejectPi", "histoRejectPi",1000, 0.51, 39, 1000, 0., 1.31 );
    SetStyleHistoTH2ForGraphs( histoRejectPi, "#it{E}_{true} (GeV/#it{c})", "#it{R}_{#pi}",
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
    histoRejectPi->GetYaxis()->SetLabelOffset(0.001);
    histoRejectPi->GetXaxis()->SetNoExponent();
    histoRejectPi->GetXaxis()->SetMoreLogLabels(kTRUE);
    TGraphErrors* graphRejectPi_90[maxcalo]   = {NULL};
    TGraphErrors* graphRejectPiSA_90[maxcalo] = {NULL};

    histoRejectPi->DrawCopy();
    TLegend* legendReject2      = GetAndSetLegend2(0.65, 0.95-textSizeLabelsRel, 0.85, 0.95-(textSizeLabelsRel +3*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,1,"",42); 
    TLegend* legendReject3      = GetAndSetLegend2(0.74, 0.95-textSizeLabelsRel, 0.94, 0.95-(textSizeLabelsRel +3*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,1,"",42); 
    TLegend* legendReject4      = GetAndSetLegend2(0.77, 0.95-textSizeLabelsRel, 0.97, 0.95-(textSizeLabelsRel +3*0.85*textSizeLabelsRel),0.85*textSizeLabelsRel,1,"",42); 
    for (Int_t icalo = 0; icalo < maxcalo; icalo++){
      graphRejectPi_90[icalo]   = new TGraphErrors(nEne,partE, pirejectM02_90[icalo]);
      while(graphRejectPi_90[icalo]->GetN()>0 && graphRejectPi_90[icalo]->GetY()[graphRejectPi_90[icalo]->GetN()-1] <= 0) graphRejectPi_90[icalo]->RemovePoint(graphRejectPi_90[icalo]->GetN()-1);
      DrawGammaSetMarkerTGraph(graphRejectPi_90[icalo], markerStyleCalo3[icalo], 1.5*markerSizeCalo[icalo], colorCalo_light2[icalo] , colorCalo_light2[icalo]);
      graphRejectPi_90[icalo]->SetLineStyle(icalo+1);
      graphRejectPi_90[icalo]->SetLineWidth(3);
      graphRejectPi_90[icalo]->Draw("samecx0");

      graphRejectPiSA_90[icalo]   = new TGraphErrors(15,partE, pirejectM20_90[icalo]);
      while(graphRejectPiSA_90[icalo]->GetN()>0 && graphRejectPiSA_90[icalo]->GetY()[graphRejectPiSA_90[icalo]->GetN()-1] <= 0) graphRejectPiSA_90[icalo]->RemovePoint(graphRejectPiSA_90[icalo]->GetN()-1);
      DrawGammaSetMarkerTGraph(graphRejectPiSA_90[icalo], markerStyleCalo3[icalo], 1.5*markerSizeCalo[icalo], colorCalo_light[icalo] , colorCalo_light[icalo]);
      graphRejectPiSA_90[icalo]->SetLineStyle(icalo+1);
      graphRejectPiSA_90[icalo]->SetLineWidth(3);
      graphRejectPiSA_90[icalo]->Draw("samecx0");

      cout << " done l:" << __LINE__ << endl;
      legendReject2->AddEntry(graphRejectPiSA_90[icalo]," ","l");
      legendReject3->AddEntry(graphRejectPi_90[icalo]," ","l");
      legendReject4->AddEntry((TObject*)0,Form("%s",caloNamePlot[icalo].Data()),"");
    }
    legendReject2->Draw();
    legendReject3->Draw();
    legendReject4->Draw();

    drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
    drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
    drawLatexAdd("#sigma_{long}^{2}",0.74,0.95-0.5*textSizeLabelsRel,textSizeLabelsRel,false,false,false);
    drawLatexAdd("#sigma_{short}^{2}",0.65,0.95-0.5*textSizeLabelsRel,textSizeLabelsRel,false,false,false);

    DrawGammaLines(0.51, 39, 1,1,2,kGray+1, 2);
  canvasReject->Update();
  canvasReject->Print(Form("%s/RejectPi_AllCalo_true.%s",outputDir.Data(),suffix.Data()));

  //********************************************************************************************************
  // create example bin for M02
  //********************************************************************************************************
  TCanvas* canvasExampleBin       = new TCanvas("canvasExampleBin", "", 200, 10, 1300, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasExampleBin,  0.105, 0.01, 0.015, 0.11);
  canvasExampleBin->SetLogy();
  cout << __LINE__ << endl;
  textSizeLabelsRel         = 58./1200;
  TH2F * histoExampleBin    = new TH2F("histoExampleBin", "histoExampleBin",1000, 0.0, 4, 1000, 0.0001, 0.805 );
  SetStyleHistoTH2ForGraphs( histoExampleBin, "#sigma_{long}^{2}", "probability",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1, 510, 505);//(#times #epsilon_{pur})
  histoExampleBin->GetYaxis()->SetLabelOffset(0.001);
  histoExampleBin->GetXaxis()->SetNoExponent();
  histoExampleBin->GetXaxis()->SetMoreLogLabels(kTRUE);
  int exampleBinCalo[maxcalo] = {3,4,5};
  Color_t colorPi = kRed+2;
  Color_t colorElec = kBlue+2;
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){

    histoExampleBin->DrawCopy();
    TLegend* legendExmplPlot           = GetAndSetLegend2(0.13, 0.885-(textSizeLabelsRel), 0.35, 0.885,0.9*textSizeLabelsRel,2,"", 42, 0.5);

    TH1D* hPionM02Plot      = (TH1D*)histM02_proj_pions[icalo][exampleBinCalo[icalo]]->Clone(Form("hPionM02Plot%s",caloName[icalo].Data()));
    TH1D* hElectronM02Plot  = (TH1D*)histM02_proj_electrons[icalo][exampleBinCalo[icalo]]->Clone(Form("hElectronM02Plot%s",caloName[icalo].Data()));
    hPionM02Plot->Scale(1/hPionM02Plot->GetEntries());
    hElectronM02Plot->Scale(1/hElectronM02Plot->GetEntries());
    DrawGammaSetMarker(hPionM02Plot, 33, 4, colorPi , colorPi);
    hPionM02Plot->Draw("samepe");
    DrawGammaSetMarker(hElectronM02Plot, 47, 3.5, colorElec , colorElec);
    hElectronM02Plot->Draw("samepe");

    legendExmplPlot->AddEntry(hPionM02Plot,"#pi^{-}","p");
    legendExmplPlot->AddEntry(hElectronM02Plot,"e^{-}","p");
    legendExmplPlot->Draw();

    drawLatexAdd(perfLabel.Data(),0.945,toplegval,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%s, w/ mat. infront",caloNamePlot[icalo].Data()),0.945,toplegval-0.05,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%1.1f< #it{#eta}< %1.1f, #it{E} = %1.1f GeV",nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1], (partE[exampleBinCalo[icalo]]+partE[exampleBinCalo[icalo]+1])/2),0.14,toplegval,textSizeLabelsRel,false,false,false);

    DrawGammaLines(M02_cut_90[icalo][exampleBinCalo[icalo]], M02_cut_90[icalo][exampleBinCalo[icalo]], 0,0.06,2,kGray+2, 2);

    canvasExampleBin->Update();
    canvasExampleBin->Print(Form("%s/M02_ExampleBin_%s.%s",outputDir.Data(),caloNamePlot[icalo].Data(),suffix.Data()));
  }

  SetStyleHistoTH2ForGraphs( histoExampleBin, "#sigma_{short}^{2}", "probability",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1, 510, 505);//(#times #epsilon_{pur})
  histoExampleBin->GetYaxis()->SetLabelOffset(0.001);
  histoExampleBin->GetXaxis()->SetNoExponent();
  histoExampleBin->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoExampleBin->GetXaxis()->SetRangeUser(0,2.5);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){

    histoExampleBin->DrawCopy();
    TLegend* legendExmplPlot           = GetAndSetLegend2(0.13, 0.885-(textSizeLabelsRel), 0.35, 0.885,0.9*textSizeLabelsRel,2,"", 42, 0.5);

    TH1D* hPionM20Plot      = (TH1D*)histM20_proj_pions[icalo][exampleBinCalo[icalo]]->Clone(Form("hPionM20Plot%s",caloName[icalo].Data()));
    TH1D* hElectronM20Plot  = (TH1D*)histM20_proj_electrons[icalo][exampleBinCalo[icalo]]->Clone(Form("hElectronM20Plot%s",caloName[icalo].Data()));
    hPionM20Plot->Scale(1/hPionM20Plot->GetEntries());
    hElectronM20Plot->Scale(1/hElectronM20Plot->GetEntries());
    DrawGammaSetMarker(hPionM20Plot, 33, 4, colorPi , colorPi);
    hPionM20Plot->Draw("samepe");
    DrawGammaSetMarker(hElectronM20Plot, 47, 3.5, colorElec , colorElec);
    hElectronM20Plot->Draw("samepe");

    legendExmplPlot->AddEntry(hPionM20Plot,"#pi^{-}","p");
    legendExmplPlot->AddEntry(hElectronM20Plot,"e^{-}","p");
    legendExmplPlot->Draw();

    drawLatexAdd(perfLabel.Data(),0.945,toplegval,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%s, w/ mat. infront",caloNamePlot[icalo].Data()),0.945,toplegval-0.05,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%1.1f< #it{#eta}< %1.1f, #it{E} = %1.1f GeV",nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1], (partE[exampleBinCalo[icalo]]+partE[exampleBinCalo[icalo]+1])/2),0.14,toplegval,textSizeLabelsRel,false,false,false);

    DrawGammaLines(M20_cut_90[icalo][exampleBinCalo[icalo]], M20_cut_90[icalo][exampleBinCalo[icalo]], 0,0.06,2,kGray+2, 2);

    canvasExampleBin->Update();
    canvasExampleBin->Print(Form("%s/M20_ExampleBin_%s.%s",outputDir.Data(),caloNamePlot[icalo].Data(),suffix.Data()));
  }
  
  // save graphs
  TFile* fileOutput = new TFile(Form("%s/output_ShowerShape.root",outputDir.Data()),"RECREATE");
  for(int icalo=0;icalo<maxcalo;icalo++){

    fileOutput->mkdir(Form("%s",caloName[icalo].Data()));
    fileOutput->cd(Form("%s",caloName[icalo].Data()));
      graphM02cutValue_90[icalo]->Write(Form("M02CutVal90_%s",caloName[icalo].Data() ));
      graphM20cutValue_90[icalo]->Write(Form("M20CutVal90_%s",caloName[icalo].Data() ));
      graphRejectPi_90[icalo]->Write(Form("RejectPiM0290_%s",caloName[icalo].Data() ));
      graphRejectPiSA_90[icalo]->Write(Form("RejectPiM2090_%s",caloName[icalo].Data() ));
    
    
  }
  // write output file   
  fileOutput->Write();
  fileOutput->Close();
}
