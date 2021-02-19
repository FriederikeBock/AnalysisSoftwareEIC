#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void resolutionBeampipe(
  Double_t jetradiusplot = 0.7,
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

  Int_t jetI = jetradiusplot==0.5 ? 0 : (jetradiusplot==0.7 ? 1 : 2);
  TString jetradius[3] = {"05","07","10"};
  Double_t etaMin[3] = {1.8,2.0,2.3};
  Double_t etaMax[3] = {3.5,3.3,3.0};
  Color_t colorRadius[3] = {kGreen+2,kOrange+2,kBlue+2};
  Marker_t markerRadius[3] = {27,28,46};

  TH2F* h_jetscale_pT[2] = {NULL};
  TH2F* h_jetscale_E[2] = {NULL};

  TTree* t_recojets_tower[2] = {NULL};
  TChain* tch_recojets_tower[2] = {NULL};

  tch_recojets_tower[0] = new TChain("ntp_recojet", "ntp_recojet");
  tch_recojets_tower[0]->Add(Form("/media/nschmidt/local/EIC_running/Fun4All_G4_FullDetector_Pythia_nobeampipe/1000/g4fwdjets_TowerFwd_%s_eval.root",jetradius[jetI].Data()));
  tch_recojets_tower[1] = new TChain("ntp_recojet", "ntp_recojet");
  tch_recojets_tower[1]->Add(Form("/media/nschmidt/local/EIC_running/Fun4All_G4_FullDetector_Pythia_beampipe/1000/g4fwdjets_TowerFwd_%s_eval.root",jetradius[jetI].Data()));
    // for(Int_t ii=0; ii<75;ii++){
    //   tch_recojets_tower[ibeam]->Add(Form("/media/nschmidt/local/EIC_afterburners/resolutionJETS/inputs/2000_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[ibeam].Data()));
    // }

  TH1D *h_JES_pT[2] = {NULL};
  TH1D *h_JER_pT[2] = {NULL};

  for(Int_t ibeam=0; ibeam<2;ibeam++){
    h_jetscale_pT[ibeam]  = new TH2F(Form("h_jetscale_pT%d",ibeam),"",40,0,40,200,-1,1);
    h_jetscale_E[ibeam]  = new TH2F(Form("h_jetscale_E%d",ibeam),"",120,0,120,200,-1,1);

    // t_recojets_tower[ibeam] = (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_afterburners/resolutionJETS/10000/g4fwdjets_TowerFwd_%s_eval.root",jetradius[ibeam].Data()), "READ"))->Get("ntp_recojet");
    if(tch_recojets_tower[ibeam]){

      Float_t rec_E,true_E;
      tch_recojets_tower[ibeam]->SetBranchAddress("e", &rec_E);
      tch_recojets_tower[ibeam]->SetBranchAddress("ge", &true_E);
      Float_t rec_pT,true_pT;
      tch_recojets_tower[ibeam]->SetBranchAddress("pt", &rec_pT);
      tch_recojets_tower[ibeam]->SetBranchAddress("gpt", &true_pT);
      Float_t rec_eta,true_eta;
      tch_recojets_tower[ibeam]->SetBranchAddress("eta", &rec_eta);
      tch_recojets_tower[ibeam]->SetBranchAddress("geta", &true_eta);

      for (Int_t i = 0; i < Int_t(tch_recojets_tower[ibeam]->GetEntries()); ++i) {
        if (tch_recojets_tower[ibeam]->LoadTree(i) < 0)
          break;

        tch_recojets_tower[ibeam]->GetEntry(i);
        if(rec_pT>1
          && (rec_eta > etaMin[ibeam] && rec_eta < etaMax[ibeam])
          ){
          h_jetscale_E[ibeam]->Fill(true_E,(rec_E-true_E)/true_E);
          h_jetscale_pT[ibeam]->Fill(true_pT,(rec_pT-true_pT)/true_pT);
        }
      }
      // h_gamma_ecalfhcalcomb_E[epart]->Scale(1 / h_gamma_ecalfhcalcomb_E[epart]->GetEntries());
    }

    h_JES_pT[ibeam] = new TH1D(Form("h_JES_pT%d",ibeam),"",40,0,40);
    h_JER_pT[ibeam] = new TH1D(Form("h_JER_pT%d",ibeam),"",40,0,40);
    for (Int_t i=1; i < h_jetscale_pT[ibeam]->GetNbinsX(); i++){
      TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[ibeam]->ProjectionY(Form("projectionYdummy%d%d",i,ibeam), i,i+1,"e");
      h_JES_pT[ibeam]->SetBinContent(i,projectionYdummy->GetMean());
      h_JES_pT[ibeam]->SetBinError(i,projectionYdummy->GetMeanError());
      h_JER_pT[ibeam]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
      h_JER_pT[ibeam]->SetBinError(i,projectionYdummy->GetRMSError()); //GetMeanError()
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
  for(Int_t ibeam=0; ibeam<2;ibeam++){
    DrawGammaSetMarker(h_JES_pT[ibeam], markerRadius[ibeam], 3, colorRadius[ibeam], colorRadius[ibeam]);
    h_JES_pT[ibeam]->Draw("same,p");
    legendJES->AddEntry(h_JES_pT[ibeam],Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot,etaMin[ibeam],etaMax[ibeam]),"p");
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
  for(Int_t ibeam=0; ibeam<2;ibeam++){
    DrawGammaSetMarker(h_JER_pT[ibeam], markerRadius[ibeam], 3, colorRadius[ibeam], colorRadius[ibeam]);
    h_JER_pT[ibeam]->Draw("same,p");
    // drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[ibeam]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[ibeam],etaMax[ibeam]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
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
  for(Int_t ibeam=0; ibeam<2;ibeam++){
    histResoDummy->Draw();

    h_jetscale_pT[ibeam]->Draw("col,same");
    DrawGammaSetMarker(h_JES_pT[ibeam], 20, 2, kMagenta+1, kMagenta+1);
    h_JES_pT[ibeam]->Draw("same");
    drawLatexAdd(collisionSystem.Data(),0.90,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[ibeam],etaMax[ibeam]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/2D_jet_pT_%s.%s", outputDir.Data(), jetradius[ibeam].Data(), suffix.Data()));
  }
  // TH2F * histResoDummyE   = new TH2F("histResoDummyE","histResoDummyE",1000,0, 119,1000,-0.5, 0.5);
  // SetStyleHistoTH2ForGraphs(histResoDummyE, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  // histResoDummyE->GetXaxis()->SetNoExponent();
  // histResoDummyE->GetYaxis()->SetNdivisions(505,kTRUE);
  // histResoDummyE->GetXaxis()->SetMoreLogLabels(kTRUE);
  // histResoDummyE->Draw();

  // h_jetscale_E[ibeam]->Draw("col,same");
  // cReso->Print(Form("%s/2D_jet_E.%s", outputDir.Data(), suffix.Data()));



}
