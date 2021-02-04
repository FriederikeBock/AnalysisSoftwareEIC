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
  // Double_t etaMin[3] = {2.5,2.5,2.5};
  // Double_t etaMax[3] = {4.0,4.0,4.0};
  Double_t etaMin[3] = {1.8,2.0,2.3};
  Double_t etaMax[3] = {3.5,3.3,3.0};
  // Double_t etaMax[3] = {3.2,3.0,2.7};
  Color_t colorRadius[3] = {kGreen+2,kOrange+2,kBlue+2};
  Marker_t markerRadius[3] = {27,28,46};

  TH2F* h_jetscale_pT[3] = {NULL};
  TH2F* h_jetscale_07_E[3] = {NULL};

  TTree* t_recojets_tower[3] = {NULL};
  TChain* tch_recojets_tower[3] = {NULL};

  for(Int_t irad=0; irad<3;irad++){

    tch_recojets_tower[irad] = new TChain("ntp_recojet", "ntp_recojet");
    tch_recojets_tower[irad]->Add(Form("/media/nschmidt/local/EIC_afterburners/resolutionJETS/10000/g4fwdjets_TowerFwd_%s_eval.root",jetradius[irad].Data()));
    for(Int_t ii=0; ii<75;ii++){
      tch_recojets_tower[irad]->Add(Form("/media/nschmidt/local/EIC_afterburners/resolutionJETS/inputs/2000_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
    }
  }

  TH1D *h_JES_pT[3] = {NULL};
  TH1D *h_JER_pT[3] = {NULL};

  for(Int_t irad=0; irad<3;irad++){
    h_jetscale_pT[irad]  = new TH2F(Form("h_jetscale_pT%d",irad),"",40,0,40,200,-1,1);
    h_jetscale_07_E[irad]  = new TH2F(Form("h_jetscale_07_E%d",irad),"",120,0,120,200,-1,1);

    // t_recojets_tower[irad] = (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_afterburners/resolutionJETS/10000/g4fwdjets_TowerFwd_%s_eval.root",jetradius[irad].Data()), "READ"))->Get("ntp_recojet");
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
        if(rec_pT>1
          && (rec_eta > etaMin[irad] && rec_eta < etaMax[irad])
          ){
          h_jetscale_07_E[irad]->Fill(true_E,(rec_E-true_E)/true_E);
          h_jetscale_pT[irad]->Fill(true_pT,(rec_pT-true_pT)/true_pT);
        }
      }
      // h_gamma_ecalfhcalcomb_E[epart]->Scale(1 / h_gamma_ecalfhcalcomb_E[epart]->GetEntries());
    }

    h_JES_pT[irad] = new TH1D(Form("h_JES_pT%d",irad),"",40,0,40);
    h_JER_pT[irad] = new TH1D(Form("h_JER_pT%d",irad),"",40,0,40);
    for (Int_t i=1; i < h_jetscale_pT[irad]->GetNbinsX(); i++){
      TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[irad]->ProjectionY(Form("projectionYdummy%d%d",i,irad), i,i+1,"e");
      h_JES_pT[irad]->SetBinContent(i,projectionYdummy->GetMean());
      h_JES_pT[irad]->SetBinError(i,projectionYdummy->GetMeanError());
      h_JER_pT[irad]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
      h_JER_pT[irad]->SetBinError(i,projectionYdummy->GetRMSError()); //GetMeanError()
    }
  }

// 2D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso, 0.1, 0.01, 0.01, 0.105);
  // cReso->SetLogz();

  TH2F* histoJESDummy   = new TH2F("histoJESDummy","histoJESDummy",1000,0, 29,1000,-0.5, 0.5);
  SetStyleHistoTH2ForGraphs(histoJESDummy, "#it{p}_{T}^{true} (GeV/#it{c})","(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histoJESDummy->GetXaxis()->SetNoExponent();
  histoJESDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histoJESDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoJESDummy->Draw();
  TLegend* legendJES  = GetAndSetLegend2(0.13, 0.83-(3*textSizeLabelsRel), 0.6, 0.83,1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
  DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
  for(Int_t irad=0; irad<3;irad++){
    DrawGammaSetMarker(h_JES_pT[irad], markerRadius[irad], 3, colorRadius[irad], colorRadius[irad]);
    h_JES_pT[irad]->Draw("same,p");
    legendJES->AddEntry(h_JES_pT[irad],Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot[irad],etaMin[irad],etaMax[irad]),"p");
  }
  legendJES->Draw();
  drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  cReso->Print(Form("%s/JES_jet_pT.%s", outputDir.Data(), suffix.Data()));


  TH2F* histoJERDummy   = new TH2F("histoJERDummy","histoJERDummy",1000,0, 29,1000,-0.01, 0.39);
  SetStyleHistoTH2ForGraphs(histoJERDummy, "#it{p}_{T}^{true} (GeV/#it{c})","#sigma(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histoJERDummy->GetXaxis()->SetNoExponent();
  histoJERDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histoJERDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoJERDummy->Draw();
  for(Int_t irad=0; irad<3;irad++){
    DrawGammaSetMarker(h_JER_pT[irad], markerRadius[irad], 3, colorRadius[irad], colorRadius[irad]);
    h_JER_pT[irad]->Draw("same,p");
    // drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  }
  legendJES->Draw();
  drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  cReso->Print(Form("%s/JER_jet_pT.%s", outputDir.Data(), suffix.Data()));


  TH2F * histResoDummy   = new TH2F("histResoDummy","histResoDummy",1000,0, 29,1000,-0.5, 0.5);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{p}_{T}^{true} (GeV/#it{c})","(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->GetZaxis()->SetRangeUser(0,500);
  for(Int_t irad=0; irad<3;irad++){
    histResoDummy->Draw();

    h_jetscale_pT[irad]->Draw("col,same");
    DrawGammaSetMarker(h_JES_pT[irad], 20, 2, kMagenta+1, kMagenta+1);
    h_JES_pT[irad]->Draw("same");
    drawLatexAdd(collisionSystem.Data(),0.90,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/2D_jet_pT_%s.%s", outputDir.Data(), jetradius[irad].Data(), suffix.Data()));
  }
  // TH2F * histResoDummyE   = new TH2F("histResoDummyE","histResoDummyE",1000,0, 119,1000,-0.5, 0.5);
  // SetStyleHistoTH2ForGraphs(histResoDummyE, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  // histResoDummyE->GetXaxis()->SetNoExponent();
  // histResoDummyE->GetYaxis()->SetNdivisions(505,kTRUE);
  // histResoDummyE->GetXaxis()->SetMoreLogLabels(kTRUE);
  // histResoDummyE->Draw();

  // h_jetscale_07_E[irad]->Draw("col,same");
  // cReso->Print(Form("%s/2D_jet_E.%s", outputDir.Data(), suffix.Data()));



}