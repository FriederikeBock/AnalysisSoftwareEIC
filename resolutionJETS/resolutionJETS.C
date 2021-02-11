#include "plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void resolutionJETS(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "e+p: 10#times250 GeV^{2}";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;

  TString jetradius[3] = {"05","07","10"};
  Double_t jetradiusplot[3] = {0.5,0.7,1.0};
  // Double_t etaMin[3] = {1.0,1.0,1.0};
  // Double_t etaMax[3] = {2.5,2.5,2.5};
  Double_t etaMin[4] = {1.5,1.5,2.0,2.5};
  Double_t etaMax[4] = {3.5,2.0,3.0,3.5};
  // Double_t etaMin[4] = {1.5,1.5,2.0,3.0};
  // Double_t etaMax[4] = {3.5,2.0,3.0,3.5};
  // Double_t etaMin[3] = {2.5,2.5,2.5};
  // Double_t etaMax[3] = {4.0,4.0,4.0};
  // Double_t etaMin[3] = {1.8,2.0,2.3};
  // Double_t etaMax[3] = {3.5,3.3,3.0};
  // Double_t etaMax[3] = {3.2,3.0,2.7};
  Color_t colorRadius[3] = {kGreen+2,kOrange+2,kBlue+2};
  Marker_t markerRadius[3] = {27,28,46};
  Color_t colorEta[4] = {kBlack,kOrange+2,kBlue+2,kGreen+2};
  Marker_t markerEta[4] = {24,33,34,47};

  TH2F* h_jetscale_pT[3][4] = {NULL};
  TH2F* h_jetscale_E[3][4] = {NULL};

  TTree* t_recojets_tower[3] = {NULL};
  TChain* tch_recojets_tower[3] = {NULL};

  for(Int_t irad=0; irad<3;irad++){

    tch_recojets_tower[irad] = new TChain("ntp_recojet", "ntp_recojet");
    tch_recojets_tower[irad]->Add(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionJETS/10000/g4fwdjets_TowerFwd_%s_eval.root",jetradius[irad].Data()));
    for(Int_t ii=0; ii<75;ii++){
      tch_recojets_tower[irad]->Add(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionJETS/inputs/2000_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
      if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()))){
        tch_recojets_tower[irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
      }
      if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()))){
        tch_recojets_tower[irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
      }
      if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()))){
        tch_recojets_tower[irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
      }
      if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()))){
        tch_recojets_tower[irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
      }
      if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()))){
        tch_recojets_tower[irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
      }
    }
    // for(Int_t ii=0; ii<26;ii++){
    //   tch_recojets_tower[irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/bad/Fun4All_G4_FullDetector_def/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
    // }
  }
Double_t xbinsDiff[16] = {0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,  20.,25.,30.,35.,40.,45.};

  TH1D *h_projectionPlot[3][4] = {NULL};
  TH1D *h_JES_pT[3][4] = {{NULL}};
  TH1D *h_JER_pT[3][4] = {{NULL}};

  for(Int_t irad=0; irad<3;irad++){
    for(Int_t etapl=0;etapl<4;etapl++){
      h_jetscale_pT[irad][etapl]  = new TH2F(Form("h_jetscale_pT%d%d",irad,etapl),"",15,xbinsDiff,200,-1,1);
      h_jetscale_E[irad][etapl]  = new TH2F(Form("h_jetscale_E%d%d",irad,etapl),"",120,0,120,200,-1,1);
    }
    // t_recojets_tower[irad] = (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionJETS/10000/g4fwdjets_TowerFwd_%s_eval.root",jetradius[irad].Data()), "READ"))->Get("ntp_recojet");
    if(tch_recojets_tower[irad]){

      Float_t rec_E,true_E;
      tch_recojets_tower[irad]->SetBranchAddress("e", &rec_E);
      tch_recojets_tower[irad]->SetBranchAddress("ge", &true_E);
      Float_t rec_pT,true_pT;
      tch_recojets_tower[irad]->SetBranchAddress("pt", &rec_pT);
      tch_recojets_tower[irad]->SetBranchAddress("gpt", &true_pT);
      Float_t rec_eta,true_eta;
      tch_recojets_tower[irad]->SetBranchAddress("eta", &rec_eta);
      tch_recojets_tower[irad]->SetBranchAddress("geta", &true_eta);

      for (Int_t i = 0; i < Int_t(tch_recojets_tower[irad]->GetEntries()); ++i) {
        if (tch_recojets_tower[irad]->LoadTree(i) < 0)
          break;

        tch_recojets_tower[irad]->GetEntry(i);
        for(Int_t etapl=0;etapl<4;etapl++){
          if(rec_pT>1
            && (rec_eta > etaMin[etapl] && rec_eta < etaMax[etapl])
            ){
            h_jetscale_E[irad][etapl]->Fill(true_E,(rec_E-true_E)/true_E);
            h_jetscale_pT[irad][etapl]->Fill(true_pT,(rec_pT-true_pT)/true_pT);
          }
        }
      }
      // h_gamma_ecalfhcalcomb_E[epart]->Scale(1 / h_gamma_ecalfhcalcomb_E[epart]->GetEntries());
    }

    for(Int_t etapl=0;etapl<4;etapl++){
      h_JES_pT[irad][etapl] = new TH1D(Form("h_JES_pT%d%d",irad,etapl),"",15,xbinsDiff);
      h_JER_pT[irad][etapl] = new TH1D(Form("h_JER_pT%d%d",irad,etapl),"",15,xbinsDiff);
      for (Int_t i=1; i < h_jetscale_pT[irad][etapl]->GetNbinsX(); i++){
        TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[irad][etapl]->ProjectionY(Form("projectionYdummy%d%d%d",i,irad,etapl), i,i+1,"e");
        // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
        h_JES_pT[irad][etapl]->SetBinContent(i,projectionYdummy->GetMean());
        h_JES_pT[irad][etapl]->SetBinError(i,projectionYdummy->GetMeanError());
        h_JER_pT[irad][etapl]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
        h_JER_pT[irad][etapl]->SetBinError(i,projectionYdummy->GetRMSError()); //GetMeanError()
        if(i==5){
          h_projectionPlot[irad][etapl] = (TH1D*)projectionYdummy->Clone(Form("h_projectionPlot_%d%d",irad,etapl));
          h_projectionPlot[irad][etapl]->Scale(1/h_projectionPlot[irad][etapl]->GetEntries());
        }
      }
    }
  }

// 2D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso, 0.11, 0.01, 0.01, 0.105);
  // cReso->SetLogz();

  TH2F* histoJESDummy   = new TH2F("histoJESDummy","histoJESDummy",1000,2.5, 34.9,1000,-0.5, 0.5);
  SetStyleHistoTH2ForGraphs(histoJESDummy, "#it{p}_{T}^{true} (GeV/#it{c})","(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histoJESDummy->GetXaxis()->SetNoExponent();
  histoJESDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histoJESDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoJESDummy->Draw();
  TLegend* legendJES  = GetAndSetLegend2(0.13, 0.83-(3*textSizeLabelsRel), 0.6, 0.83,1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
  DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
  for(Int_t irad=0; irad<3;irad++){
    DrawGammaSetMarker(h_JES_pT[irad][0], markerRadius[irad], 3, colorRadius[irad], colorRadius[irad]);
    h_JES_pT[irad][0]->Draw("same,p");
    legendJES->AddEntry(h_JES_pT[irad][0],Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot[irad],etaMin[0],etaMax[0]),"p");
  }
  legendJES->Draw();
  drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  cReso->Print(Form("%s/JES_jet_pT.%s", outputDir.Data(), suffix.Data()));


  TH2F* histoJERDummy   = new TH2F("histoJERDummy","histoJERDummy",1000,2.5, 34.9,1000,-0.01, 0.39);
  SetStyleHistoTH2ForGraphs(histoJERDummy, "#it{p}_{T}^{true} (GeV/#it{c})","#sigma(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.99);
  histoJERDummy->GetXaxis()->SetNoExponent();
  histoJERDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histoJERDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoJERDummy->Draw();
  for(Int_t irad=0; irad<3;irad++){
    DrawGammaSetMarker(h_JER_pT[irad][0], markerRadius[irad], 3, colorRadius[irad], colorRadius[irad]);
    h_JER_pT[irad][0]->Draw("same,p");
    // drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  }
  legendJES->Draw();
  DrawGammaLines(2.5,34.9,0.1,0.1,2,kGray+1,7);
  drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  cReso->Print(Form("%s/JER_jet_pT.%s", outputDir.Data(), suffix.Data()));

  for(Int_t irad=0; irad<3;irad++){
    histoJERDummy->Draw();
    TLegend* legendJES2  = GetAndSetLegend2(0.13, 0.83-(4*textSizeLabelsRel), 0.6, 0.83,1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    for(Int_t etapl=0; etapl<4;etapl++){
      // if(etapl==3){h_JER_pT[irad][etapl]->Rebin(2);h_JER_pT[irad][etapl]->Scale(0.5);}
      DrawGammaSetMarker(h_JER_pT[irad][etapl], markerEta[etapl], 3, colorEta[etapl], colorEta[etapl]);
      h_JER_pT[irad][etapl]->Draw("same,p");
      legendJES2->AddEntry(h_JER_pT[irad][etapl],Form("%1.1f < #eta_{jet} < %1.1f",etaMin[etapl],etaMax[etapl]),"p");
      // drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      // drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }
    legendJES2->Draw();
    DrawGammaLines(2.5,34.9,0.1,0.1,2,kGray+1,7);
    drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} < %1.1f, HCAL+ECAL jets",jetradiusplot[irad]),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cReso->Print(Form("%s/JER_jet_pT_R%1.1f.%s", outputDir.Data(),jetradiusplot[irad], suffix.Data()));
  }

  for(Int_t etapl=0; etapl<4;etapl++){
    histoJERDummy->Draw();
    TLegend* legendJES2  = GetAndSetLegend2(0.13, 0.83-(3*textSizeLabelsRel), 0.6, 0.83,1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    for(Int_t irad=0; irad<3;irad++){
      DrawGammaSetMarker(h_JER_pT[irad][etapl], markerRadius[irad], 3, colorRadius[irad], colorRadius[irad]);
      h_JER_pT[irad][etapl]->Draw("same,p");
      legendJES2->AddEntry(h_JER_pT[irad][etapl],Form("#it{R} < %1.1f",jetradiusplot[irad]),"p");
      // drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      // drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }
    legendJES2->Draw();
    DrawGammaLines(2.5,34.9,0.1,0.1,2,kGray+1,7);
    drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, %1.1f < #eta_{jet} < %1.1f, HCAL+ECAL jets",etaMin[etapl],etaMax[etapl]),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cReso->Print(Form("%s/JER_jet_pT_Eta%d.%s", outputDir.Data(),etapl, suffix.Data()));
  }

  TH2F * histResoDummy   = new TH2F("histResoDummy","histResoDummy",1000,5, 34.9,1000,-0.5, 0.5);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{p}_{T}^{true} (GeV/#it{c})","(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->GetZaxis()->SetRangeUser(0,500);
  for(Int_t irad=0; irad<3;irad++){
    histResoDummy->Draw();

    h_jetscale_pT[irad][0]->Draw("col,same");
    DrawGammaSetMarker(h_JES_pT[irad][0], 20, 2, kMagenta+1, kMagenta+1);
    h_JES_pT[irad][0]->Draw("same");
    drawLatexAdd(collisionSystem.Data(),0.90,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[0],etaMax[0]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    histResoDummy->Draw("same,axis");
    cReso->Print(Form("%s/2D_jet_pT_%s.%s", outputDir.Data(), jetradius[irad].Data(), suffix.Data()));
  }
  // TH2F * histResoDummyE   = new TH2F("histResoDummyE","histResoDummyE",1000,0, 119,1000,-0.5, 0.5);
  // SetStyleHistoTH2ForGraphs(histResoDummyE, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  // histResoDummyE->GetXaxis()->SetNoExponent();
  // histResoDummyE->GetYaxis()->SetNdivisions(505,kTRUE);
  // histResoDummyE->GetXaxis()->SetMoreLogLabels(kTRUE);
  // histResoDummyE->Draw();

  // h_jetscale_E[irad]->Draw("col,same");
  // cReso->Print(Form("%s/2D_jet_E.%s", outputDir.Data(), suffix.Data()));

// 2D PLOT
  TCanvas* cSingleSlice = new TCanvas("cSingleSlice","",0,0,1100,1000);
  DrawGammaCanvasSettings( cSingleSlice, 0.11, 0.01, 0.002, 0.105);
  // cSingleSlice->SetLogz();

  TH2F* histoJESSliceDummy   = new TH2F("histoJESSliceDummy","histoJESSliceDummy",1000,-0.69, 0.09,1000,0., 0.079);
  SetStyleHistoTH2ForGraphs(histoJESSliceDummy, "(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}","norm. counts ", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  histoJESSliceDummy->GetYaxis()->SetNoExponent();
  histoJESSliceDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  // histoJESSliceDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  for(Int_t irad=0; irad<3;irad++){

    histoJESSliceDummy->Draw();
    // TLegend* legendJES  = GetAndSetLegend2(0.13, 0.83-(3*textSizeLabelsRel), 0.6, 0.83,1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    // DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
      TLegend* legendJES3  = GetAndSetLegend2(0.13, 0.80-(5*textSizeLabelsRel), 0.6, 0.80-(1*textSizeLabelsRel),1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    for(Int_t etapl=0; etapl<4;etapl++){
      DrawGammaSetMarker(h_projectionPlot[irad][etapl], markerEta[etapl], 3, colorEta[etapl], colorEta[etapl]);
      h_projectionPlot[irad][etapl]->Draw("same,p");
      // legendJES->AddEntry(h_JES_pT[irad][0],Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot[irad],etaMin[0],etaMax[0]),"p");
      legendJES3->AddEntry(h_projectionPlot[irad][etapl],Form("%1.1f < #eta_{jet} < %1.1f",etaMin[etapl],etaMax[etapl]),"p");
    }
    legendJES3->Draw();
    drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f, HCAL+ECAL (tower) jets",jetradiusplot[irad]),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("8.0 < #it{p}_{T}^{rec} < 10 GeV/#it{c}"),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cSingleSlice->Print(Form("%s/JES_Slice_Plot_R%d.%s", outputDir.Data(),irad, suffix.Data()));
  }

  for(Int_t etapl=0; etapl<4;etapl++){
    histoJESSliceDummy->Draw();
    // TLegend* legendJES  = GetAndSetLegend2(0.13, 0.83-(3*textSizeLabelsRel), 0.6, 0.83,1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    // DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
      TLegend* legendJES3  = GetAndSetLegend2(0.13, 0.80-(4*textSizeLabelsRel), 0.6, 0.80-(1*textSizeLabelsRel),1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    for(Int_t irad=0; irad<3;irad++){
      DrawGammaSetMarker(h_projectionPlot[irad][etapl], markerRadius[irad], 3, colorRadius[irad], colorRadius[irad]);
      h_projectionPlot[irad][etapl]->Draw("same,p");
      // legendJES->AddEntry(h_JES_pT[irad][0],Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot[irad],etaMin[0],etaMax[0]),"p");
      legendJES3->AddEntry(h_projectionPlot[irad][etapl],Form("#it{R} = %1.1f",jetradiusplot[irad]),"p");
    }
    legendJES3->Draw();
    drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL (tower) jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[etapl],etaMax[etapl]),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("8.0 < #it{p}_{T,jet}^{rec} < 10 GeV/#it{c}"),0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cSingleSlice->Print(Form("%s/JES_Slice_Plot2_Eta%d.%s", outputDir.Data(),etapl, suffix.Data()));
  }



}