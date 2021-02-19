#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))


void runbasictreeloop(TChain* inpChain, Int_t iInp);

void studies_leak_acceptance_pipe(
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
  Color_t colorInput[4] = {kGreen+2,kOrange+2,kBlue+2,kRed+2};
  Marker_t markerInput[4] = {27,28,46,20};

  TString inputFolderNames[4] = {"Fun4All_G4_FullDetector_Pythia_nobeampipe","Fun4All_G4_FullDetector_Pythia_beampipe","Fun4All_G4_HighAccDetector_Pythia","Fun4All_G4_FullEtaLeakTestDetector_Pythia"};

  TChain* tch_jets_07[4]  = {NULL};
  TH2F* h_jetscale_pT[4]  = {NULL};
  TH1D *h_JES_pT[4]       = {NULL};
  TH1D *h_JER_pT[4]       = {NULL};
  Double_t etaMin = 2.3;
  Double_t etaMax = 4.0;

  TChain* tch_cluster_hcal[4]        = {NULL};
  TChain* tch_cluster_ecal[4]        = {NULL};
  TChain* tch_shower_hcal[4]         = {NULL};
  TChain* tch_shower_ecal[4]         = {NULL};
  const Int_t nEne = 25;
  const static Double_t partE[]   = {1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0,  10., 11., 13., 15., 20., 25., 30., 40., 50., 60., 70., 80., 150.};
  TH1F* h_pion_fhcal_E[nEne] 	    = {NULL};
  TH1F* h_pion_fhcal_E[nEne] 	    = {NULL};

  for(Int_t iInp=0; iInp<4;iInp++){
    tch_cluster_hcal[iInp] = new TChain("ntp_cluster", "ntp_cluster");
    tch_cluster_hcal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/%s/1000/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data()));

    tch_cluster_ecal[iInp] = new TChain("ntp_cluster", "ntp_cluster");
    tch_cluster_ecal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/%s/1000/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data()));

    tch_shower_hcal[iInp] = new TChain("ntp_gshower", "ntp_gshower");
    tch_shower_hcal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/%s/1000/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data()));

    tch_shower_ecal[iInp] = new TChain("ntp_gshower", "ntp_gshower");
    tch_shower_ecal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/%s/1000/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data()));

    for(Int_t ii=0; ii<100;ii++){
      if(!gSystem->AccessPathName(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii))){
        tch_cluster_hcal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii));
        tch_shower_hcal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii));
      }
      if(!gSystem->AccessPathName(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),ii))){
        tch_cluster_hcal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),ii));
        tch_shower_hcal[iInp]->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),ii));
      }
    }

    // JETS
    tch_jets_07[iInp] = new TChain("ntp_recojet", "ntp_recojet");
    tch_jets_07[iInp]->Add(Form("/media/nschmidt/local/EIC_running/%s/1000/g4fwdjets_TowerFwd_07_eval.root",inputFolderNames[iInp].Data()));
    for(Int_t ii=0; ii<100;ii++){
      if(!gSystem->AccessPathName(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/g4fwdjets_TowerFwd_07_eval.root",inputFolderNames[iInp].Data(),ii))){
        tch_jets_07[iInp]->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/500_%d/g4fwdjets_TowerFwd_07_eval.root",inputFolderNames[iInp].Data(),ii));
      }
    }

    h_jetscale_pT[iInp]  = new TH2F(Form("h_jetscale_pT%d",iInp),"",40,0,40,200,-1,1);

    if(tch_jets_07[iInp]){
      Float_t rec_pT,true_pT;
      tch_jets_07[iInp]->SetBranchAddress("pt", &rec_pT);
      tch_jets_07[iInp]->SetBranchAddress("gpt", &true_pT);
      Float_t rec_eta,true_eta;
      tch_jets_07[iInp]->SetBranchAddress("eta", &rec_eta);
      tch_jets_07[iInp]->SetBranchAddress("geta", &true_eta);

      for (Int_t i = 0; i < Int_t(tch_jets_07[iInp]->GetEntries()); ++i) {
        if (tch_jets_07[iInp]->LoadTree(i) < 0)
          break;

        tch_jets_07[iInp]->GetEntry(i);
        if(rec_pT>1
          && (rec_eta > etaMin && rec_eta < etaMax)
          ){
          h_jetscale_pT[iInp]->Fill(true_pT,(rec_pT-true_pT)/true_pT);
        }
      }
    }

    h_JES_pT[iInp] = new TH1D(Form("h_JES_pT%d",iInp),"",40,0,40);
    h_JER_pT[iInp] = new TH1D(Form("h_JER_pT%d",iInp),"",40,0,40);
    for (Int_t i=1; i < h_jetscale_pT[iInp]->GetNbinsX(); i++){
      TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[iInp]->ProjectionY(Form("projectionYdummy%d%d",i,iInp), i,i+1,"e");
      h_JES_pT[iInp]->SetBinContent(i,projectionYdummy->GetMean());
      h_JES_pT[iInp]->SetBinError(i,projectionYdummy->GetMeanError());
      h_JER_pT[iInp]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
      h_JER_pT[iInp]->SetBinError(i,projectionYdummy->GetRMSError()); //GetMeanError()
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
  for(Int_t irad=0; irad<4;irad++){
    DrawGammaSetMarker(h_JES_pT[irad], markerInput[irad], 3, colorInput[irad], colorInput[irad]);
    h_JES_pT[irad]->Draw("same,p");
    legendJES->AddEntry(h_JES_pT[irad],Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot,etaMin,etaMax),"p");
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
  for(Int_t irad=0; irad<4;irad++){
    DrawGammaSetMarker(h_JER_pT[irad], markerInput[irad], 3, colorInput[irad], colorInput[irad]);
    h_JER_pT[irad]->Draw("same,p");
    // drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  }
  legendJES->Draw();
  drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  cReso->Print(Form("%s/JER_jet_pT.%s", outputDir.Data(), suffix.Data()));


}
