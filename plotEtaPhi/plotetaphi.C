#include "plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);

void plotetaphi(
    Int_t plotEvent = 7,
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;

  Int_t granularity = 4;
  TH2F* h_IEtaIPhiMapEvt_HCAL 	    = NULL;
  TH2F* h_IEtaIPhiMapEvt_HCAL_Ecut 	    = NULL;
  TH2F* h_EtaPhiMapEvt_HCAL 	    = NULL;
  TChain* t_towers_HCAL = new TChain("ntp_tower", "ntp_tower");
  // t_towers_HCAL->Add("/media/nschmidt/local/EIC_running/cluster_output/Fun4All_G4_FullDetector_Pythia_beampipe/plotting/JET_EIC_FST1/G4EICDetector_g4fhcal_eval.root");
  t_towers_HCAL->Add("/media/nschmidt/local/EIC_running/Fun4All_G4_FullDetector_MaxHighGranHCAL_Pions_beampipe/3000_50.0/G4EICDetector_g4fhcal_eval.root");
  // t_towers_HCAL->Add("/media/nschmidt/local/EIC_running/Fun4All_G4_FullDetector_HighGranHCAL_Pions_beampipe/3000_50.0/G4EICDetector_g4fhcal_eval.root");
  // t_towers_HCAL->Add("/media/nschmidt/local/EIC_running/Fun4All_G4_FullDetector_DefaultGranHCAL_beampipe/3000_50.0/G4EICDetector_g4fhcal_eval.root");
  if(t_towers_HCAL){
    h_IEtaIPhiMapEvt_HCAL 	= new TH2F("h_IEtaIPhiMapEvt_HCAL", "",  granularity*55, -0.5, granularity*55-0.5, granularity*55, -0.5, granularity*55-0.5);
    h_IEtaIPhiMapEvt_HCAL_Ecut 	= new TH2F("h_IEtaIPhiMapEvt_HCAL_Ecut", "",  granularity*55, -0.5, granularity*55-0.5, granularity*55, -0.5, granularity*55-0.5);
    h_EtaPhiMapEvt_HCAL 	= new TH2F("h_EtaPhiMapEvt_HCAL", "", granularity*55, -3.14, 3.14,  granularity*55,1.2, 4.0);
    Float_t phi_HCAL,eta_HCAL,iphi_HCAL,ieta_HCAL,e_HCAL,event;
    t_towers_HCAL->SetBranchAddress("event", &event);
    t_towers_HCAL->SetBranchAddress("phi", &phi_HCAL);
    t_towers_HCAL->SetBranchAddress("eta", &eta_HCAL);
    t_towers_HCAL->SetBranchAddress("iphi", &iphi_HCAL);
    t_towers_HCAL->SetBranchAddress("ieta", &ieta_HCAL);
    t_towers_HCAL->SetBranchAddress("e", &e_HCAL);
    Int_t nentries = Int_t(t_towers_HCAL->GetEntries());
    Float_t lastEvt = 0;
    for (Int_t i = 0; i < nentries; ++i) {
      if (t_towers_HCAL->LoadTree(i) < 0)
        break;
      t_towers_HCAL->GetEntry(i);
      if(event!=plotEvent) continue;

      h_IEtaIPhiMapEvt_HCAL->Fill(iphi_HCAL,ieta_HCAL,e_HCAL);
      if(e_HCAL>0.1)h_IEtaIPhiMapEvt_HCAL_Ecut->Fill(iphi_HCAL,ieta_HCAL,e_HCAL);

      h_EtaPhiMapEvt_HCAL->Fill(phi_HCAL,eta_HCAL,e_HCAL);
    }
  }

  TH2F* h_IEtaIPhiMapEvt_ECAL 	    = NULL;
  TH2F* h_EtaPhiMapEvt_ECAL 	    = NULL;
  TChain* t_towers_ECAL = new TChain("ntp_tower", "ntp_tower");
  t_towers_ECAL->Add("/media/nschmidt/local/EIC_running/cluster_output/Fun4All_G4_FullDetector_Pythia_beampipe/plotting/JET_EIC_FST1/G4EICDetector_g4femc_eval.root");
  if(t_towers_ECAL){
    h_IEtaIPhiMapEvt_ECAL 	= new TH2F("h_IEtaIPhiMapEvt_ECAL", "",  66, -0.5, 65.5, 66, -0.5, 65.5);
    Float_t phi_ECAL,eta_ECAL,iphi_ECAL,ieta_ECAL,e_ECAL;
    t_towers_ECAL->SetBranchAddress("phi", &phi_ECAL);
    t_towers_ECAL->SetBranchAddress("eta", &eta_ECAL);
    t_towers_ECAL->SetBranchAddress("iphi", &iphi_ECAL);
    t_towers_ECAL->SetBranchAddress("ieta", &ieta_ECAL);
    t_towers_ECAL->SetBranchAddress("e", &e_ECAL);
    Int_t nentries = Int_t(t_towers_ECAL->GetEntries());
    Float_t lastEvt = 0;
    for (Int_t i = 0; i < nentries; ++i) {
      if (t_towers_ECAL->LoadTree(i) < 0)
        break;
      t_towers_ECAL->GetEntry(i);
      h_IEtaIPhiMapEvt_ECAL->Fill(iphi_ECAL,ieta_ECAL,e_ECAL);
    }
  }



  // 2D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1000,950);
  DrawGammaCanvasSettings( cReso, 0.1, 0.03, 0.01, 0.105);
  // cReso->SetLogz();

  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_HCAL, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.2,1.4);
  h_IEtaIPhiMapEvt_HCAL->GetZaxis()->SetTitle("tower #it{E} (GeV)");
  h_IEtaIPhiMapEvt_HCAL->GetZaxis()->SetTitleOffset(1.2);
  h_IEtaIPhiMapEvt_HCAL->Draw("lego2");
    drawLatexAdd("FHCAL, e-p: 10#times250 GeV^{2}",0.97,0.02,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // cReso->Print(Form("%s/IEtaIPhiHCAL.%s", outputDir.Data(), suffix.Data()));

  // 1D PLOT
  cReso->SetLeftMargin(0.12);
  cReso->SetRightMargin(0.1);
  cReso->SetBottomMargin(0.1);
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_HCAL, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  h_IEtaIPhiMapEvt_HCAL->Draw("colz");
   TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
   el4->SetFillStyle(0);
   el4->SetLineColor(kGray+2);
   el4->SetLineWidth(2);
   el4->SetLineStyle(7);
   el4->Draw("same");
   if(granularity==1){
  DrawGammaLines(25.5, 28.5, 28.5, 28.5, 2, kGray+3, 7);
  DrawGammaLines(25.5, 28.5, 25.5, 25.5, 2, kGray+3, 7);
  DrawGammaLines(25.5, 25.5, 25.5, 28.5, 2, kGray+3, 7);
  DrawGammaLines(28.5, 28.5, 25.5, 28.5, 2, kGray+3, 7);
   }
   if(granularity==2){
  DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*28.5-1, granularity*28.5-1, 2, kGray+3, 7);
  DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*25.5-1, granularity*25.5-1, 2, kGray+3, 7);
  DrawGammaLines(granularity*25.5-1, granularity*25.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
  DrawGammaLines(granularity*28.5-1, granularity*28.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
   }
  if(granularity==4){
  TEllipse *el5 = new TEllipse(granularity*27.0-2,granularity*27.0-2,granularity*1.5,granularity*1.5,0,360,0);
   el5->SetFillStyle(0);
   el5->SetLineColor(kGray+3);
   el5->SetLineWidth(2);
   el5->SetLineStyle(7);
   el5->Draw("same");
   }
    drawLatexAdd("FHCAL, e-p: 10#times250 GeV^{2}",0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/IEtaIPhiHCAL_1D%d.%s", outputDir.Data(),plotEvent, suffix.Data()));

  // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_HCAL_Ecut, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  h_IEtaIPhiMapEvt_HCAL_Ecut->Draw("colz");
   el4->Draw("same");
   if(granularity==1){
  DrawGammaLines(25.5, 28.5, 28.5, 28.5, 2, kGray+3, 7);
  DrawGammaLines(25.5, 28.5, 25.5, 25.5, 2, kGray+3, 7);
  DrawGammaLines(25.5, 25.5, 25.5, 28.5, 2, kGray+3, 7);
  DrawGammaLines(28.5, 28.5, 25.5, 28.5, 2, kGray+3, 7);
   }
   if(granularity==2){
  DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*28.5-1, granularity*28.5-1, 2, kGray+3, 7);
  DrawGammaLines(granularity*25.5-1, granularity*28.5-1, granularity*25.5-1, granularity*25.5-1, 2, kGray+3, 7);
  DrawGammaLines(granularity*25.5-1, granularity*25.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
  DrawGammaLines(granularity*28.5-1, granularity*28.5-1, granularity*25.5-1, granularity*28.5-1, 2, kGray+3, 7);
   }
   if(granularity==4){
  TEllipse *el5 = new TEllipse(granularity*27.0-2,granularity*27.0-2,granularity*1.5,granularity*1.5,0,360,0);
   el5->SetFillStyle(0);
   el5->SetLineColor(kGray+3);
   el5->SetLineWidth(2);
   el5->SetLineStyle(7);
   el5->Draw("same");
   }
    drawLatexAdd("FHCAL, e-p: 10#times250 GeV^{2}",0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/IEtaIPhiHCAL_1D%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));

  SetStyleHistoTH2ForGraphs(h_EtaPhiMapEvt_HCAL, "tower #phi","tower #eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  h_EtaPhiMapEvt_HCAL->GetZaxis()->SetTitle("tower #it{E}");
  h_EtaPhiMapEvt_HCAL->Draw("lego2");
  // cReso->Print(Form("%s/EtaPhiHCAL.%s", outputDir.Data(), suffix.Data()));

  // 1D PLOT
  h_EtaPhiMapEvt_HCAL->Draw("col");
  // cReso->Print(Form("%s/EtaPhiHCAL_1D.%s", outputDir.Data(), suffix.Data()));


  // 2D PLOT
  cReso = new TCanvas("cReso","",0,0,1000,1000);
  DrawGammaCanvasSettings( cReso, 0.1, 0.03, 0.01, 0.105);
  // cReso->SetLogz();

  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_ECAL, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.2,1.4);
  h_IEtaIPhiMapEvt_ECAL->GetZaxis()->SetTitle("tower #it{E} (GeV)");
  h_IEtaIPhiMapEvt_ECAL->GetZaxis()->SetTitleOffset(1.2);
  h_IEtaIPhiMapEvt_ECAL->Draw("lego2");
    drawLatexAdd("FECAL, e-p: 10#times250 GeV^{2}",0.97,0.02,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // cReso->Print(Form("%s/IEtaIPhiECAL.%s", outputDir.Data(), suffix.Data()));

  // 1D PLOT
  h_IEtaIPhiMapEvt_ECAL->Draw("col");
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
    drawLatexAdd("FECAL, e-p: 10#times250 GeV^{2}",0.93,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // cReso->Print(Form("%s/EtaPhiECAL_1D.%s", outputDir.Data(), suffix.Data()));


}


void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs)
{
  if(numInputs<5){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 2, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(2, 2);
  }else if(numInputs<7){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(3, 2);
  }else if(numInputs<10){
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
