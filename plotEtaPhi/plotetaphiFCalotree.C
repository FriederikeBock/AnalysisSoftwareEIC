#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);

void plotetaphiFCalotree(
    TString nameFile  = "",
    Int_t plotEvent   = 1,
    Double_t minE     = 40.,
    TString suffix            = "pdf"
){
  int towersx = 55;
  int towersxFEMC = 64;
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
  TH2F* h_IEtaIPhiMapEvt_FHCAL 	    = NULL;
  TH2F* h_IEtaIPhiMapEvt_FHCAL_Ecut 	    = NULL;
  TH2F* h_IEtaIPhiMapEvt_FEMC 	    = NULL;
  TH2F* h_IEtaIPhiMapEvt_FEMC_Ecut 	    = NULL;
  TChain* t_towers_FHCAL = new TChain("event_tree", "event_tree");
  t_towers_FHCAL->Add(nameFile.Data());
  const int _maxNTowersFHCAL  = 50 * 50;
  const int _maxNTowersFEMC   = 150 * 150;
  Double_t totEnergy      = 0.;
  Double_t totEnergyFEMC  = 0.;
  Bool_t foundIt          = kFALSE;
  if(t_towers_FHCAL){
    h_IEtaIPhiMapEvt_FHCAL 	= new TH2F("h_IEtaIPhiMapEvt_FHCAL", "",  towersx+12, -6.5, towersx+5.5,  towersx+12, -6.5, towersx+5.5);
    h_IEtaIPhiMapEvt_FHCAL_Ecut 	= new TH2F("h_IEtaIPhiMapEvt_FHCAL_Ecut", "",  towersx+12, -6.5, towersx+5.5, towersx+12, -6.5, towersx+5.5);
    h_IEtaIPhiMapEvt_FEMC 	= new TH2F("h_IEtaIPhiMapEvt_FEMC", "",  towersxFEMC+12, -6.5, towersxFEMC+5.5,  towersxFEMC+12, -6.5, towersxFEMC+5.5);
    h_IEtaIPhiMapEvt_FEMC_Ecut 	= new TH2F("h_IEtaIPhiMapEvt_FEMC_Ecut", "",  towersxFEMC+12, -6.5, towersxFEMC+5.5, towersxFEMC+12, -6.5, towersxFEMC+5.5);
    int _nTowers_FHCAL;
    float* _tower_FHCAL_E            = new float[_maxNTowersFHCAL];
    int* _tower_FHCAL_iEta         = new int[_maxNTowersFHCAL];
    int* _tower_FHCAL_iPhi         = new int[_maxNTowersFHCAL];
    int* _tower_FHCAL_trueID       = new int[_maxNTowersFHCAL];
    int _nTowers_FEMC;
    float* _tower_FEMC_E            = new float[_maxNTowersFEMC];
    int* _tower_FEMC_iEta         = new int[_maxNTowersFEMC];
    int* _tower_FEMC_iPhi         = new int[_maxNTowersFEMC];
    int* _tower_FEMC_trueID       = new int[_maxNTowersFEMC];
    t_towers_FHCAL->SetBranchAddress("tower_FHCAL_N",                &_nTowers_FHCAL);
    t_towers_FHCAL->SetBranchAddress("tower_FHCAL_E",                _tower_FHCAL_E);
    t_towers_FHCAL->SetBranchAddress("tower_FHCAL_iEta",             _tower_FHCAL_iEta);
    t_towers_FHCAL->SetBranchAddress("tower_FHCAL_iPhi",             _tower_FHCAL_iPhi);
    t_towers_FHCAL->SetBranchAddress("tower_FEMC_N",                &_nTowers_FEMC);
    t_towers_FHCAL->SetBranchAddress("tower_FEMC_E",                _tower_FEMC_E);
    t_towers_FHCAL->SetBranchAddress("tower_FEMC_iEta",             _tower_FEMC_iEta);
    t_towers_FHCAL->SetBranchAddress("tower_FEMC_iPhi",             _tower_FEMC_iPhi);
    Int_t nentries = Int_t(t_towers_FHCAL->GetEntries());
    for (Int_t i = 0; i < nentries && !foundIt; i++) {
      if (t_towers_FHCAL->LoadTree(i) < 0)
        break;
      if(i!=plotEvent) continue;
      totEnergy     = 0.;
      totEnergyFEMC = 0.;
      t_towers_FHCAL->GetEntry(i);
      for(int itwr=0;itwr<_nTowers_FHCAL;itwr++){
        h_IEtaIPhiMapEvt_FHCAL->Fill(_tower_FHCAL_iEta[itwr],_tower_FHCAL_iPhi[itwr],_tower_FHCAL_E[itwr]*2);
        if(_tower_FHCAL_E[itwr]>0.1)h_IEtaIPhiMapEvt_FHCAL_Ecut->Fill(_tower_FHCAL_iEta[itwr],_tower_FHCAL_iPhi[itwr],_tower_FHCAL_E[itwr]*2);
        totEnergy=totEnergy+_tower_FHCAL_E[itwr]*2;
        cout << itwr << "\t" << _tower_FHCAL_E[itwr]*2 << endl;
      }
      if (totEnergy < minE){
        cout << i << "\t FHCAL: " << totEnergy  << endl;
      }
      for(int itwr=0;itwr<_nTowers_FEMC;itwr++){
        h_IEtaIPhiMapEvt_FEMC->Fill(_tower_FEMC_iEta[itwr],_tower_FEMC_iPhi[itwr],_tower_FEMC_E[itwr]);
        if(_tower_FEMC_E[itwr]>0.1)h_IEtaIPhiMapEvt_FEMC_Ecut->Fill(_tower_FEMC_iEta[itwr],_tower_FEMC_iPhi[itwr],_tower_FEMC_E[itwr]);
        totEnergyFEMC=totEnergyFEMC+_tower_FEMC_E[itwr];
        cout << itwr << "\t" << _tower_FEMC_E[itwr]*2 << endl;
      }
      if (totEnergy < minE){
        h_IEtaIPhiMapEvt_FHCAL->Clear();
        h_IEtaIPhiMapEvt_FHCAL_Ecut->Clear();
        h_IEtaIPhiMapEvt_FEMC->Clear();
        h_IEtaIPhiMapEvt_FEMC_Ecut->Clear();
        cout << i <<  "\t FEMC: "<< totEnergyFEMC << endl;
        plotEvent++;
      } else {
        foundIt = kTRUE;
      }
    }
  }
  cout << "\t FHCAL: " << totEnergy << "\t FEMC: " << totEnergyFEMC << endl;
  
  // 2D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1000,1000);
  DrawGammaCanvasSettings( cReso, 0.1, 0.015, 0.015, 0.1);
  // cReso->SetLogz();

  // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_FHCAL,"tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  h_IEtaIPhiMapEvt_FHCAL->Draw("col");
   TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
   el4->SetFillStyle(0);
   el4->SetLineColor(kGray+2);
   el4->SetLineWidth(2);
   el4->SetLineStyle(7);
   el4->Draw("same");
   if(granularity==1){
      DrawGammaLines(24.5, 29.5, 29.5, 29.5, 2, kGray+3, 7);
      DrawGammaLines(24.5, 29.5, 24.5, 24.5, 2, kGray+3, 7);
      DrawGammaLines(24.5, 24.5, 24.5, 29.5, 2, kGray+3, 7);
      DrawGammaLines(29.5, 29.5, 24.5, 29.5, 2, kGray+3, 7);
   }
  //  if(granularity==2){
  // DrawGammaLines(granularity*24.5-1, granularity*29.5-1, granularity*29.5-1, granularity*29.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*24.5-1, granularity*29.5-1, granularity*24.5-1, granularity*24.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*24.5-1, granularity*24.5-1, granularity*24.5-1, granularity*29.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*29.5-1, granularity*29.5-1, granularity*24.5-1, granularity*29.5-1, 2, kGray+3, 7);
  //  }
  // if(granularity==4){
  // TEllipse *el5 = new TEllipse(granularity*27.0-2,granularity*27.0-2,granularity*1.5,granularity*1.5,0,360,0);
  //  el5->SetFillStyle(0);
  //  el5->SetLineColor(kGray+3);
  //  el5->SetLineWidth(2);
  //  el5->SetLineStyle(7);
  //  el5->Draw("same");
  //  }
    drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  cReso->Print(Form("%s/IEtaIPhiFHCAL_1D%d.%s", outputDir.Data(),plotEvent, suffix.Data()));
  cReso->Print(Form("%s/IEtaIPhiFHCAL_1D.%s", outputDir.Data(), suffix.Data()));

  // // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_FHCAL_Ecut, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad,  0.9,0.9);
  h_IEtaIPhiMapEvt_FHCAL_Ecut->Draw("col");
   el4->Draw("same");
   if(granularity==1){
      DrawGammaLines(24.5, 29.5, 29.5, 29.5, 2, kGray+3, 7);
      DrawGammaLines(24.5, 29.5, 24.5, 24.5, 2, kGray+3, 7);
      DrawGammaLines(24.5, 24.5, 24.5, 29.5, 2, kGray+3, 7);
      DrawGammaLines(29.5, 29.5, 24.5, 29.5, 2, kGray+3, 7);
   }
//    if(granularity==2){
  // DrawGammaLines(granularity*24.5-1, granularity*29.5-1, granularity*29.5-1, granularity*29.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*24.5-1, granularity*29.5-1, granularity*24.5-1, granularity*24.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*24.5-1, granularity*24.5-1, granularity*24.5-1, granularity*29.5-1, 2, kGray+3, 7);
  // DrawGammaLines(granularity*29.5-1, granularity*29.5-1, granularity*24.5-1, granularity*29.5-1, 2, kGray+3, 7);
  //  }
  //  if(granularity==4){
  // TEllipse *el5 = new TEllipse(granularity*27.0-2,granularity*27.0-2,granularity*1.5,granularity*1.5,0,360,0);
  //  el5->SetFillStyle(0);
  //  el5->SetLineColor(kGray+3);
  //  el5->SetLineWidth(2);
  //  el5->SetLineStyle(7);
  //  el5->Draw("same");
  //  }
    drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  cReso->Print(Form("%s/IEtaIPhiFHCAL_1D_Ecut.%s", outputDir.Data(), suffix.Data()));
  cReso->Print(Form("%s/IEtaIPhiFHCAL_1D%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));

  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_FEMC, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad,  0.9,0.9);
  h_IEtaIPhiMapEvt_FEMC->Draw("col");
   el4 = new TEllipse(32.0,32.0,32,32,0,360,0);
   el4->SetFillStyle(0);
   el4->SetLineColor(kGray+2);
   el4->SetLineWidth(2);
   el4->SetLineStyle(7);
   el4->Draw("same");
    DrawGammaLines(28.5, 35.5, 35.5, 35.5, 2, kGray+2, 7);
    DrawGammaLines(28.5, 35.5, 28.5, 28.5, 2, kGray+2, 7);
    DrawGammaLines(28.5, 28.5, 28.5, 35.5, 2, kGray+2, 7);
    DrawGammaLines(35.5, 35.5, 28.5, 35.5, 2, kGray+2, 7);
    drawLatexAdd("FEMC, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#it{E}^{FEMC}_{tot} = %.1f GeV", totEnergyFEMC),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  cReso->Print(Form("%s/IEtaIPhiFEMC_1D%d.%s", outputDir.Data(),plotEvent, suffix.Data()));
  cReso->Print(Form("%s/IEtaIPhiFEMC_1D.%s", outputDir.Data(), suffix.Data()));

  // // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_FEMC_Ecut, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  h_IEtaIPhiMapEvt_FEMC_Ecut->Draw("col");
   el4->Draw("same");
   if(granularity==1){
      DrawGammaLines(24.5, 29.5, 29.5, 29.5, 2, kGray+3, 7);
      DrawGammaLines(24.5, 29.5, 24.5, 24.5, 2, kGray+3, 7);
      DrawGammaLines(24.5, 24.5, 24.5, 29.5, 2, kGray+3, 7);
      DrawGammaLines(29.5, 29.5, 24.5, 29.5, 2, kGray+3, 7);
   }
    drawLatexAdd("FEMC, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#it{E}^{FEMC}_{tot} = %.1f GeV", totEnergyFEMC),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  cReso->Print(Form("%s/IEtaIPhiFEMC_1D_Ecut.%s", outputDir.Data(), suffix.Data()));
  cReso->Print(Form("%s/IEtaIPhiFEMC_1D%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));
    
}

