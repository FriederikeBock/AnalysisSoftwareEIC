#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);

void plotetaphitree(
    Int_t plotEvent = 7,
    TString suffix            = "pdf"
){
  int towersx = 440/0.15;
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;

  Int_t granularity = 1;
  TH1F* h_IDvsScint_DRCALO 	    = NULL;
  TH1F* h_IDvsCerenkov_DRCALO 	    = NULL;
  TH2F* h_IEtaIPhiMapEvt_DRCALO 	    = NULL;
  TH2F* h_IEtaIPhiMapEvt_DRCALO_Ecut 	    = NULL;
  TChain* t_towers_DRCALO = new TChain("event_tree", "event_tree");
  t_towers_DRCALO->Add("/media/nschmidt/local/EIC_running/DualReadout/DualReadout_DRCALOONLY_DRLARGE_SimplePion_testingHighRes/eventtree.root");
  const int _maxNTowersDR = towersx * towersx;
  if(t_towers_DRCALO){
    h_IDvsScint_DRCALO 	= new TH1F("h_IDvsScint_DRCALO", "",  500, 1250, 1750);
    h_IDvsCerenkov_DRCALO 	= new TH1F("h_IDvsCerenkov_DRCALO", "",  500, 1250, 1750);
    h_IEtaIPhiMapEvt_DRCALO 	= new TH2F("h_IEtaIPhiMapEvt_DRCALO", "",  towersx, -0.5, towersx-0.5, towersx, -0.5, towersx-0.5);
    h_IEtaIPhiMapEvt_DRCALO_Ecut 	= new TH2F("h_IEtaIPhiMapEvt_DRCALO_Ecut", "",  towersx, -0.5, towersx-0.5, towersx, -0.5, towersx-0.5);
    int _nTowers_DRCALO;
   float* _tower_DRCALO_E            = new float[_maxNTowersDR];
    int* _tower_DRCALO_iEta         = new int[_maxNTowersDR];
    int* _tower_DRCALO_iPhi         = new int[_maxNTowersDR];
    int* _tower_DRCALO_trueID       = new int[_maxNTowersDR];
    int* _tower_DRCALO_NCerenkov       = new int[_maxNTowersDR];
    int* _tower_DRCALO_NScint       = new int[_maxNTowersDR];
    t_towers_DRCALO->SetBranchAddress("tower_DRCALO_N",                &_nTowers_DRCALO);
    t_towers_DRCALO->SetBranchAddress("tower_DRCALO_E",                _tower_DRCALO_E);
    t_towers_DRCALO->SetBranchAddress("tower_DRCALO_iEta",             _tower_DRCALO_iEta);
    t_towers_DRCALO->SetBranchAddress("tower_DRCALO_iPhi",             _tower_DRCALO_iPhi);
    t_towers_DRCALO->SetBranchAddress("tower_DRCALO_NCerenkov",             _tower_DRCALO_NCerenkov);
    t_towers_DRCALO->SetBranchAddress("tower_DRCALO_NScint",             _tower_DRCALO_NScint);
    Int_t nentries = Int_t(t_towers_DRCALO->GetEntries());
    for (Int_t i = 0; i < nentries; ++i) {
      if (t_towers_DRCALO->LoadTree(i) < 0)
        break;
      if(i!=plotEvent) continue;
      t_towers_DRCALO->GetEntry(i);
      for(int itwr=0;itwr<_nTowers_DRCALO;itwr++){
        if(_tower_DRCALO_iPhi[itwr]==1200 && _tower_DRCALO_NCerenkov[itwr]>0){
          cout << _tower_DRCALO_iEta[itwr] << "\t" << _tower_DRCALO_NCerenkov[itwr] << endl;
          h_IDvsCerenkov_DRCALO->Fill(_tower_DRCALO_iEta[itwr]);
        }
        if(_tower_DRCALO_iPhi[itwr]==1200 && _tower_DRCALO_NScint[itwr]>0)h_IDvsScint_DRCALO->Fill(_tower_DRCALO_iEta[itwr]);
        h_IEtaIPhiMapEvt_DRCALO->Fill(_tower_DRCALO_iEta[itwr],_tower_DRCALO_iPhi[itwr],_tower_DRCALO_E[itwr]);
        if(_tower_DRCALO_E[itwr]>0.1)h_IEtaIPhiMapEvt_DRCALO_Ecut->Fill(_tower_DRCALO_iEta[itwr],_tower_DRCALO_iPhi[itwr],_tower_DRCALO_E[itwr]);
      }
    }
  }



  // 2D PLOT
  TCanvas* cResoiEta = new TCanvas("cResoiEta","",0,0,4000,1000);
  DrawGammaCanvasSettings( cResoiEta, 0.1, 0.03, 0.01, 0.105);
  cResoiEta->SetLogy();
  // 1D PLOT
  cResoiEta->SetLeftMargin(0.12);
  cResoiEta->SetRightMargin(0.1);
  cResoiEta->SetBottomMargin(0.1);
  SetStyleHistoTH1ForGraphs(h_IDvsCerenkov_DRCALO,"tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  h_IDvsCerenkov_DRCALO->Draw("col");

    drawLatexAdd("Dual-Readout, e-p: 10#times250 GeV^{2}",0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cResoiEta->Print(Form("%s/iEtavsCherenkov%d.%s", outputDir.Data(),plotEvent, suffix.Data()));
  cResoiEta->Print(Form("%s/iEtavsCherenkov_1D.%s", outputDir.Data(), suffix.Data()));

  SetStyleHistoTH1ForGraphs(h_IDvsScint_DRCALO,"tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  h_IDvsScint_DRCALO->Draw("col");

    drawLatexAdd("Dual-Readout, e-p: 10#times250 GeV^{2}",0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cResoiEta->Print(Form("%s/iEtavsScint%d.%s", outputDir.Data(),plotEvent, suffix.Data()));
  cResoiEta->Print(Form("%s/iEtavsScint_1D.%s", outputDir.Data(), suffix.Data()));




  // 2D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1000,950);
  DrawGammaCanvasSettings( cReso, 0.1, 0.03, 0.01, 0.105);
  // cReso->SetLogz();

  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_DRCALO, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.2,1.4);
  h_IEtaIPhiMapEvt_DRCALO->GetZaxis()->SetTitle("tower #it{E} (GeV)");
  h_IEtaIPhiMapEvt_DRCALO->GetZaxis()->SetTitleOffset(1.2);
  h_IEtaIPhiMapEvt_DRCALO->Draw("lego2");
    drawLatexAdd("DRCALO, e-p: 10#times250 GeV^{2}",0.97,0.02,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // cReso->Print(Form("%s/IEtaIPhiDRCALO.%s", outputDir.Data(), suffix.Data()));

  // 1D PLOT
  cReso->SetLeftMargin(0.12);
  cReso->SetRightMargin(0.1);
  cReso->SetBottomMargin(0.1);
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_DRCALO,"tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  h_IEtaIPhiMapEvt_DRCALO->Draw("colz");
  //  TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
  //  el4->SetFillStyle(0);
  //  el4->SetLineColor(kGray+2);
  //  el4->SetLineWidth(2);
  //  el4->SetLineStyle(7);
  //  el4->Draw("same");
  //  if(granularity==1){
  // DrawGammaLines(25.5, 28.5, 28.5, 28.5, 2, kGray+3, 7);
  // DrawGammaLines(25.5, 28.5, 25.5, 25.5, 2, kGray+3, 7);
  // DrawGammaLines(25.5, 25.5, 25.5, 28.5, 2, kGray+3, 7);
  // DrawGammaLines(28.5, 28.5, 25.5, 28.5, 2, kGray+3, 7);
  //  }
  //  if(granularity==2){
  // DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*28.5-1, granularity*28.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*25.5-1, granularity*25.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*25.5-1, granularity*25.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*28.5-1, granularity*28.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
  //  }
  // if(granularity==4){
  // TEllipse *el5 = new TEllipse(granularity*27.0-2,granularity*27.0-2,granularity*1.5,granularity*1.5,0,360,0);
  //  el5->SetFillStyle(0);
  //  el5->SetLineColor(kGray+3);
  //  el5->SetLineWidth(2);
  //  el5->SetLineStyle(7);
  //  el5->Draw("same");
  //  }
    drawLatexAdd("Dual-Readout, e-p: 10#times250 GeV^{2}",0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/IEtaIPhiDRCALO_1D%d.%s", outputDir.Data(),plotEvent, suffix.Data()));
  cReso->Print(Form("%s/IEtaIPhiDRCALO_1D.%s", outputDir.Data(), suffix.Data()));

  // // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_DRCALO_Ecut, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  h_IEtaIPhiMapEvt_DRCALO_Ecut->Draw("colz");
  //  el4->Draw("same");
  //  if(granularity==1){
  // DrawGammaLines(25.5, 28.5, 28.5, 28.5, 2, kGray+3, 7);
  // DrawGammaLines(25.5, 28.5, 25.5, 25.5, 2, kGray+3, 7);
  // DrawGammaLines(25.5, 25.5, 25.5, 28.5, 2, kGray+3, 7);
  // DrawGammaLines(28.5, 28.5, 25.5, 28.5, 2, kGray+3, 7);
  //  }
  //  if(granularity==2){
  // DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*28.5-1, granularity*28.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*25.5-1, granularity*25.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*25.5-1, granularity*25.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*28.5-1, granularity*28.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
  //  }
  //  if(granularity==4){
  // TEllipse *el5 = new TEllipse(granularity*27.0-2,granularity*27.0-2,granularity*1.5,granularity*1.5,0,360,0);
  //  el5->SetFillStyle(0);
  //  el5->SetLineColor(kGray+3);
  //  el5->SetLineWidth(2);
  //  el5->SetLineStyle(7);
  //  el5->Draw("same");
  //  }
    drawLatexAdd("DRCALO, e-p: 10#times250 GeV^{2}",0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/IEtaIPhiDRCALO_1D_Ecut.%s", outputDir.Data(), suffix.Data()));
  cReso->Print(Form("%s/IEtaIPhiDRCALO_1D%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));

  // SetStyleHistoTH2ForGraphs(h_EtaPhiMapEvt_DRCALO, "tower #phi","tower #eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  // h_EtaPhiMapEvt_DRCALO->GetZaxis()->SetTitle("tower #it{E}");
  // h_EtaPhiMapEvt_DRCALO->Draw("lego2");
  // // cReso->Print(Form("%s/EtaPhiDRCALO.%s", outputDir.Data(), suffix.Data()));

  // // 1D PLOT
  // h_EtaPhiMapEvt_DRCALO->Draw("col");
  // // cReso->Print(Form("%s/EtaPhiDRCALO_1D.%s", outputDir.Data(), suffix.Data()));


  // // 2D PLOT
  // cReso = new TCanvas("cReso","",0,0,1000,1000);
  // DrawGammaCanvasSettings( cReso, 0.1, 0.03, 0.01, 0.105);
  // // cReso->SetLogz();

  // SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_ECAL, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.2,1.4);
  // h_IEtaIPhiMapEvt_ECAL->GetZaxis()->SetTitle("tower #it{E} (GeV)");
  // h_IEtaIPhiMapEvt_ECAL->GetZaxis()->SetTitleOffset(1.2);
  // h_IEtaIPhiMapEvt_ECAL->Draw("lego2");
  //   drawLatexAdd("FECAL, e-p: 10#times250 GeV^{2}",0.97,0.02,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // // cReso->Print(Form("%s/IEtaIPhiECAL.%s", outputDir.Data(), suffix.Data()));

  // // 1D PLOT
  // h_IEtaIPhiMapEvt_ECAL->Draw("col");
  //  el4 = new TEllipse(32.0,32.0,32,32,0,360,0);
  //  el4->SetFillStyle(0);
  //  el4->SetLineColor(kGray+2);
  //  el4->SetLineWidth(2);
  //  el4->SetLineStyle(7);
  //  el4->Draw("same");
  // DrawGammaLines(28.5, 35.5, 35.5, 35.5, 2, kGray+2, 7);
  // DrawGammaLines(28.5, 35.5, 28.5, 28.5, 2, kGray+2, 7);
  // DrawGammaLines(28.5, 28.5, 28.5, 35.5, 2, kGray+2, 7);
  // DrawGammaLines(35.5, 35.5, 28.5, 35.5, 2, kGray+2, 7);
  //   drawLatexAdd("FECAL, e-p: 10#times250 GeV^{2}",0.93,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // // cReso->Print(Form("%s/EtaPhiECAL_1D.%s", outputDir.Data(), suffix.Data()));


}

