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
void eoverpplotting_paper(
    TString inputFileName     = "central_allculsters.root",
    TString suffix            = "pdf"
){  
  TString clusterizerName = "MA clusters";
  TString collisionsSys     = "PythiaMB";
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
  TString dateForOutput                       = ReturnDateStringForOutput();

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
  
  TString outputDir                 = Form("plotsEoverP_paper/%s",dateForOutput.Data());
  // if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFile                      = new TFile(inputFileName.Data());
  if(inputFile->IsZombie()){
    cout << "inputFile not found!" << endl;
  }

  //************************** Read data **************************************************
//   const int nEne = 15;
//   const static Double_t partE[]   = {   0.7, 1.0, 1.5, 2, 3, 4, 6, 10,15,20,25,  30,35,  40, 50};
  const int nEne = 28;
  const static Double_t partE[]   = {   0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2, 3,
                                        4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                        14, 15, 16, 17, 18, 19, 20, 25,  30};
  Double_t partEcenter[nEne];
  for (int i = 0; i < nEne; i++){
    partEcenter[i]  = (partE[i]+partE[i+1])/2.;
  }
  
  int exampleBinCalo[maxcalo] = {3,4,5};
  Color_t colorPi = kRed+2;
  Color_t colorElec = kBlue+2;
  double toplegval = 0.90;

  TH1D* hPionEoPplot[maxcalo] = {nullptr};
  TH1D* hElectronEoPplot[maxcalo] = {nullptr};
  TH1D* hPionEoPplot_rec[maxcalo] = {nullptr};
  TH1D* hElectronEoPplot_rec[maxcalo] = {nullptr};
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    hPionEoPplot[icalo]    = (TH1D*)inputFile->Get(Form("hPionEoPplot_%s_%d",caloName[icalo].Data(),exampleBinCalo[icalo]));
    hElectronEoPplot[icalo]    = (TH1D*)inputFile->Get(Form("hElectronEoPplot_%s_%d",caloName[icalo].Data(),exampleBinCalo[icalo]));
    hPionEoPplot_rec[icalo]    = (TH1D*)inputFile->Get(Form("hPionEoPplot_rec_%s_%d",caloName[icalo].Data(),exampleBinCalo[icalo]));
    hElectronEoPplot_rec[icalo]    = (TH1D*)inputFile->Get(Form("hElectronEoPplot_rec_%s_%d",caloName[icalo].Data(),exampleBinCalo[icalo]));
  }

  TCanvas* canvasExampleBin       = new TCanvas("canvasExampleBin", "", 200, 10, 1300, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasExampleBin,  0.105, 0.01, 0.015, 0.09);
  // canvasExampleBin->SetLogy(1);
  // canvasExampleBin->SetLogx(1);
  cout << __LINE__ << endl;
  textSizeLabelsRel        = 58./1200;
  TH2F * histoExampleBin;
  histoExampleBin                = new TH2F("histoExampleBin", "histoExampleBin",1000, 0.0, 1.23, 1000, -0.002, 0.23 );
  SetStyleHistoTH2ForGraphs( histoExampleBin, "#it{E}/#it{p}", "probability",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1, 510, 505);//(#times #epsilon_{pur})
  histoExampleBin->GetYaxis()->SetLabelOffset(0.001);
  histoExampleBin->GetXaxis()->SetNoExponent();
  histoExampleBin->GetXaxis()->SetMoreLogLabels(true);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){

    histoExampleBin->DrawCopy();
    TLegend* legendExmplPlot           = GetAndSetLegend2(0.13, 0.885-(2*textSizeLabelsRel), 0.55, 0.885,0.9*textSizeLabelsPixel,2);
    // TLegend* legendPi0Merge2           = GetAndSetLegend2(0.75, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    // legendPi0Merge3           = GetAndSetLegend2(0.70, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    int ieta = nEta;

    // histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][exampleBinCalo[icalo]] 
    // histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][exampleBinCalo[icalo]]
    DrawGammaSetMarker(hPionEoPplot[icalo], 33, 4, colorPi , colorPi);
    hPionEoPplot[icalo]->Draw("samepe");
    DrawGammaSetMarker(hElectronEoPplot[icalo], 47, 3.5, colorElec , colorElec);
    hElectronEoPplot[icalo]->Draw("samepe");

    DrawGammaSetMarker(hPionEoPplot_rec[icalo], 27, 4, colorPi-6 , colorPi-6);
    hPionEoPplot_rec[icalo]->Draw("samepe");
    DrawGammaSetMarker(hElectronEoPplot_rec[icalo], 46, 3.5, colorElec-6 , colorElec-6);
    hElectronEoPplot_rec[icalo]->Draw("samepe");

    legendExmplPlot->AddEntry(hPionEoPplot[icalo],"#pi^{-} (true #it{p})","p");
    legendExmplPlot->AddEntry(hElectronEoPplot[icalo],"e^{-} (true #it{p})","p");
    legendExmplPlot->AddEntry(hPionEoPplot_rec[icalo],"#pi^{-} (rec. #it{p})","p");
    legendExmplPlot->AddEntry(hElectronEoPplot_rec[icalo],"e^{-} (rec. #it{p})","p");

    legendExmplPlot->Draw();
    drawLatexAdd(perfLabel.Data(),0.945,toplegval,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%s, w/ mat. infront",caloNamePlot[icalo].Data()),0.945,toplegval-0.05,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%1.1f< #it{#eta}< %1.1f, #it{E} = %1.1f GeV",nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1], (partE[exampleBinCalo[icalo]]+partE[exampleBinCalo[icalo]+1])/2),0.14,toplegval,textSizeLabelsRel,false,false,false);

    // DrawGammaLines(EP_cut[1][imatch][icalo][ieta][exampleBinCalo[icalo]], EP_cut[1][imatch][icalo][nEta][exampleBinCalo[icalo]], 0,0.1,2,kGray+2, 2);

    canvasExampleBin->Update();
    canvasExampleBin->Print(Form("%s/EoP_ExampleBin_%s.%s",outputDir.Data(),caloNamePlot[icalo].Data(),suffix.Data()));
  }



    Double_t arrayBoundariesTMEffEtaX1_4[4];
    Double_t arrayBoundariesTMEffEtaY1_4[3];
    Double_t relativeMarginsTMEffEtaX[3];
    Double_t relativeMarginsTMEffEtaY[3];
    textSizeLabelsPixel             = 750*0.05;
    ReturnCorrectValuesForCanvasScaling(2250,750, 3, 1,0.04, 0.003, 0.01,0.095,arrayBoundariesTMEffEtaX1_4,arrayBoundariesTMEffEtaY1_4,relativeMarginsTMEffEtaX,relativeMarginsTMEffEtaY);

    TCanvas* canvasTMEffEtaPaperPlot     = new TCanvas("canvasTMEffEtaPaperPlot","",0,0,2250,750);  // gives the page size
    DrawGammaCanvasSettings( canvasTMEffEtaPaperPlot,  0.05, 0.01, 0.01,0.095);
    // canvasTMEffEtaPaperPlot->SetLogy(1);

    TPad* padTMEffEtaEEMC               = new TPad("padTMEffEtaEEMC", "", arrayBoundariesTMEffEtaX1_4[0], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaEEMC, relativeMarginsTMEffEtaX[0], relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaY[0], relativeMarginsTMEffEtaY[2]);
    padTMEffEtaEEMC->Draw();

    TPad* padTMEffEtaBEMC                = new TPad("padTMEffEtaBEMC", "", arrayBoundariesTMEffEtaX1_4[1], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaBEMC, relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaY[0], relativeMarginsTMEffEtaY[2]);
    padTMEffEtaBEMC->Draw();

    TPad* padTMEffEtaFEMC                = new TPad("padTMEffEtaFEMC", "", arrayBoundariesTMEffEtaX1_4[2], arrayBoundariesTMEffEtaY1_4[2], arrayBoundariesTMEffEtaX1_4[3], arrayBoundariesTMEffEtaY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padTMEffEtaFEMC, relativeMarginsTMEffEtaX[1], relativeMarginsTMEffEtaX[2], relativeMarginsTMEffEtaY[0], relativeMarginsTMEffEtaY[2]);
    padTMEffEtaFEMC->Draw();

    // padTMEffEtaEEMC->cd();
    // padTMEffEtaEEMC->SetLogy();

    Double_t margin                 = relativeMarginsTMEffEtaX[0]*2.7*1350;
    double textsizeLabelsTMEffEta    = 0;
    double textsizeFacTMEffEta       = 0;
    if (padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) < padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1())){
        textsizeLabelsTMEffEta         = (Double_t)textSizeLabelsPixel/padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) ;
        textsizeFacTMEffEta            = (Double_t)1./padTMEffEtaEEMC->XtoPixel(padTMEffEtaEEMC->GetX2()) ;
    } else {
        textsizeLabelsTMEffEta         = (Double_t)textSizeLabelsPixel/padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1());
        textsizeFacTMEffEta            = (Double_t)1./padTMEffEtaEEMC->YtoPixel(padTMEffEtaEEMC->GetY1());
    }
    histoExampleBin                = new TH2F("histoExampleBin", "histoExampleBin",1000, 0.001, 1.25, 1000, -0.002, 0.23 );
    SetStyleHistoTH2ForGraphs( histoExampleBin, "#it{E}/#it{p}", "probability",
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1, 510, 505);//(#times #epsilon_{pur})
    histoExampleBin->GetYaxis()->SetLabelOffset(0.001);
    histoExampleBin->GetXaxis()->SetNoExponent();
    histoExampleBin->GetXaxis()->SetMoreLogLabels(true);
    TLegend* legendExmplPlot           = GetAndSetLegend2(0.05, 0.95-(2*textSizeLabelsRel), 0.55, 0.95,0.9*textSizeLabelsPixel,2);
    for (Int_t icalo = 0; icalo < 3; icalo++){
      if(icalo==1){
        padTMEffEtaBEMC->cd();
      } else if(icalo==2){
        padTMEffEtaFEMC->cd();
      } else {
        padTMEffEtaEEMC->cd();
      }

      histoExampleBin->DrawCopy();
      if(icalo==0){
        double startleg = 0.91;
        drawLatexAdd(perfLabel,0.16,startleg,textSizeLabelsRel,false,false,false);
        drawLatexAdd(Form("single e^{#pm}"),0.16, startleg-(1)*textSizeLabelsRel*1.1,textSizeLabelsRel,false,false,false);

        drawLatexAdd( caloNamePlot[icalo].Data(),0.95,0.91,textSizeLabelsRel,false,false,true);
        drawLatexAdd(Form("%1.1f< #it{#eta}< %1.1f",nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1]),0.95,0.86,textSizeLabelsRel,false,false,true);
        drawLatexAdd(Form("#it{E} = %1.1f GeV",(partE[exampleBinCalo[icalo]]+partE[exampleBinCalo[icalo]+1])/2),0.95,0.81,textSizeLabelsRel,false,false,true);
        drawLatexAdd( "a)",0.16,0.15,textSizeLabelsRel,false,false,false);
      } else if(icalo==1){
        drawLatexAdd( caloNamePlot[icalo].Data(),0.95,0.91,textSizeLabelsRel,false,false,true);
        drawLatexAdd(Form("%1.1f< #it{#eta}< %1.1f",nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1]),0.95,0.86,textSizeLabelsRel,false,false,true);
        drawLatexAdd(Form("#it{E} = %1.1f GeV",(partE[exampleBinCalo[icalo]]+partE[exampleBinCalo[icalo]+1])/2),0.95,0.81,textSizeLabelsRel,false,false,true);
        drawLatexAdd( "b)",0.05,0.15,textSizeLabelsRel,false,false,false);
      } else {
        drawLatexAdd( caloNamePlot[icalo].Data(),0.93,0.91,textSizeLabelsRel,false,false,true);
        drawLatexAdd(Form("%1.1f< #it{#eta}< %1.1f",nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1]),0.93,0.86,textSizeLabelsRel,false,false,true);
        drawLatexAdd(Form("#it{E} = %1.1f GeV",(partE[exampleBinCalo[icalo]]+partE[exampleBinCalo[icalo]+1])/2),0.93,0.81,textSizeLabelsRel,false,false,true);
        drawLatexAdd("c)",0.05,0.15,textSizeLabelsRel,false,false,false);
      }

      DrawGammaSetMarker(hPionEoPplot[icalo], 33, 2, colorPi , colorPi);
      hPionEoPplot[icalo]->Draw("samepe");
      DrawGammaSetMarker(hElectronEoPplot[icalo], 47, 1.8, colorElec , colorElec);
      hElectronEoPplot[icalo]->Draw("samepe");

      DrawGammaSetMarker(hPionEoPplot_rec[icalo], 27, 2, colorPi-6 , colorPi-6);
      hPionEoPplot_rec[icalo]->Draw("samepe");
      DrawGammaSetMarker(hElectronEoPplot_rec[icalo], 46, 1.8, colorElec-6 , colorElec-6);
      hElectronEoPplot_rec[icalo]->Draw("samepe");

      if(icalo==1){
        legendExmplPlot->AddEntry(hPionEoPplot[icalo],"#pi^{-} (true #it{p})","p");
        legendExmplPlot->AddEntry(hElectronEoPplot[icalo],"e^{-} (true #it{p})","p");
        legendExmplPlot->AddEntry(hPionEoPplot_rec[icalo],"#pi^{-} (rec. #it{p})","p");
        legendExmplPlot->AddEntry(hElectronEoPplot_rec[icalo],"e^{-} (rec. #it{p})","p");
        legendExmplPlot->Draw();
      }
      // DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      // Int_t icurrPID = 1;
      // Int_t nActiveEta            = maxNEtaBinsFull[region[icalo]]+1;
      // legendEffiE[icalo]      = GetAndSetLegend2(0.4, 0.94-(nActiveEta/2*0.85*textSizeLabelsRel), 0.95, 0.94,0.85*textSizeLabelsRel, 2, "", 42, 0.2);
      // for (Int_t iEta=minEtaBinCaloDis[region[icalo]]; iEta<maxEtaBinCaloDis[region[icalo]]+1;iEta++){
      //   int iEtaPlot = iEta==maxEtaBinCaloDis[region[icalo]] ? nEta : iEta;
      //   DrawGammaSetMarker(h_TMeffi_recSE_MCE[icalo][icurrPID][iEta], markerStyleEta[iEtaPlot], markerSizeEta[iEtaPlot], colorEta[iEtaPlot], colorEta[iEtaPlot]);
      //   h_TMeffi_recSE_MCE[icalo][icurrPID][iEta]->Draw("same");

      //   if (iEta == maxEtaBinCaloDis[region[icalo]] )
      //     legendEffiE[icalo]->AddEntry(h_TMeffi_recSE_MCE[icalo][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[minEtaBinCaloDis[region[icalo]]],partEta[maxEtaBinCaloDis[region[icalo]]]),"p");
      //   else 
      //     legendEffiE[icalo]->AddEntry(h_TMeffi_recSE_MCE[icalo][icurrPID][iEta],Form("%1.1f < #it{#eta} < %1.1f",partEta[iEta],partEta[iEta+1]),"p");
      //   legendEffiE[icalo]->Draw();
      // }
      histoExampleBin->Draw("same,axis");
    }
    canvasTMEffEtaPaperPlot->SaveAs(Form("%s/EoP_ExampleBin_PaperPlot.%s",outputDir.Data(), suffix.Data()));



}
