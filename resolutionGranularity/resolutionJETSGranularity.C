#include "plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void resolutionJETSGranularity(
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

  Color_t colorInput[4] = {kOrange+2,kMagenta+2,kCyan+2,kRed+2};
  Marker_t markerInput[4] = {28,24,46,20};

  TString jetradius[3] = {"05","07","10"};
  Double_t jetradiusplot[3] = {0.5,0.7,1.0};
  // Double_t etaMin[3] = {1.0,1.0,1.0};
  // Double_t etaMin[3] = {2.0,2.0,2.0};
  Double_t etaMin[3] = {2.5,2.5,2.5};
  // Double_t etaMin[3] = {3.0,3.0,3.0};
  // Double_t etaMax[3] = {4.0,4.0,4.0};
  Double_t etaMax[3] = {3.5,3.5,3.5};
  // Double_t etaMin[3] = {1.8,2.0,2.3};
  // Double_t etaMax[3] = {3.5,3.3,3.0};
  // Double_t etaMax[3] = {3.2,3.0,2.7};
  Color_t colorRadius[3] = {kGreen+2,kOrange+2,kBlue+2};
  Marker_t markerRadius[3] = {27,28,46};

  TH2F* h_jetscale_pT[3][3] = {NULL};
  TH2F* h_jetscale_E[3][3] = {NULL};
  TH2F* h_jet_dEtadPhi[3][3] = {NULL};

  TTree* t_recojets_tower[3][3] = {NULL};
  TChain* tch_recojets_tower[3][3] = {NULL};

  TString inputFolderNames[3] = {"Fun4All_G4_FullDetector_DefaultGranHCAL_PythiaJET_beampipe","Fun4All_G4_FullDetector_HighGranHCAL_PythiaJET_beampipe","Fun4All_G4_FullDetector_MaxHighGranHCAL_PythiaJET_beampipe"};
  TString namePlotInput[3] = {"default","2x granularity","4x granularity"};
Double_t xbinsDiff[16] = {0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,  20.,25.,30.,35.,40.,45.};

  TH1D *h_JES_pT[3][3] = {NULL};
  TH1D *h_JER_pT[3][3] = {NULL};
  for(Int_t iInp=0; iInp<3;iInp++){
    for(Int_t irad=0; irad<3;irad++){

      tch_recojets_tower[iInp][irad] = new TChain("ntp_recojet", "ntp_recojet");
      for(Int_t ii=0; ii<26;ii++){
        if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/500_%d/g4fwdjets_TowerFwd_%s_eval.root",inputFolderNames[iInp].Data(),ii,jetradius[irad].Data()))){
        tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/500_%d/g4fwdjets_TowerFwd_%s_eval.root",inputFolderNames[iInp].Data(),ii,jetradius[irad].Data()));
        }
        if(iInp==0){
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
        }
        if(iInp==1){
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC2x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC2x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC2x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC2x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC2x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC2x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC2x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC2x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC2x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC2x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC2x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC2x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
        }
        if(iInp==2){
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC4x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC4x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC4x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG10-HC4x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC4x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC4x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC4x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG20-HC4x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC4x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC4x/2000_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
          if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC4x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data())))
            tch_recojets_tower[iInp][irad]->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/Fun4All_G4_FullDetector_JTTRG30-HC4x/300_-1_%d/g4fwdjets_TowerFwd_%s_eval.root",ii,jetradius[irad].Data()));
        }
      }
    }

    for(Int_t irad=0; irad<3;irad++){
      h_jetscale_pT[iInp][irad]  = new TH2F(Form("h_jetscale_pT%d",irad),"",15,xbinsDiff,200,-1,1);
      h_jetscale_E[iInp][irad]  = new TH2F(Form("h_jetscale_E%d",irad),"",120,0,120,200,-1,1);
      h_jet_dEtadPhi[iInp][irad]  = new TH2F(Form("h_jet_dEtadPhi%d",irad),"",200,-0.5,0.5,200,-0.5,0.5);

      if(tch_recojets_tower[iInp][irad]){

      Float_t rec_E,true_E;
      tch_recojets_tower[iInp][irad]->SetBranchAddress("e", &rec_E);
      tch_recojets_tower[iInp][irad]->SetBranchAddress("ge", &true_E);
      Float_t rec_pT,true_pT;
      tch_recojets_tower[iInp][irad]->SetBranchAddress("pt", &rec_pT);
      tch_recojets_tower[iInp][irad]->SetBranchAddress("gpt", &true_pT);
      Float_t rec_eta,true_eta;
      tch_recojets_tower[iInp][irad]->SetBranchAddress("eta", &rec_eta);
      tch_recojets_tower[iInp][irad]->SetBranchAddress("geta", &true_eta);
      Float_t rec_phi,true_phi;
      tch_recojets_tower[iInp][irad]->SetBranchAddress("phi", &rec_phi);
      tch_recojets_tower[iInp][irad]->SetBranchAddress("gphi", &true_phi);

      for (Int_t i = 0; i < Int_t(tch_recojets_tower[iInp][irad]->GetEntries()); ++i) {
        if (tch_recojets_tower[iInp][irad]->LoadTree(i) < 0)
          break;

        tch_recojets_tower[iInp][irad]->GetEntry(i);
        if(rec_pT>1
          && (rec_eta > etaMin[irad] && rec_eta < etaMax[irad])
          ){
          h_jetscale_E[iInp][irad]->Fill(true_E,(rec_E-true_E)/true_E);
          h_jetscale_pT[iInp][irad]->Fill(true_pT,(rec_pT-true_pT)/true_pT);
          h_jet_dEtadPhi[iInp][irad]->Fill(rec_eta-true_eta,rec_phi-true_phi);
        }
      }
      h_jet_dEtadPhi[iInp][irad]->Scale(1 / h_jet_dEtadPhi[iInp][irad]->GetEntries());
    }

    h_JES_pT[iInp][irad] = new TH1D(Form("h_JES_pT%d",irad),"",15,xbinsDiff);
    h_JER_pT[iInp][irad] = new TH1D(Form("h_JER_pT%d",irad),"",15,xbinsDiff);
    for (Int_t i=1; i < h_jetscale_pT[iInp][irad]->GetNbinsX(); i++){
      TH1D* projectionYdummy = (TH1D*)h_jetscale_pT[iInp][irad]->ProjectionY(Form("projectionYdummy%d%d",i,irad), i,i+1,"e");
      h_JES_pT[iInp][irad]->SetBinContent(i,projectionYdummy->GetMean());
      h_JES_pT[iInp][irad]->SetBinError(i,projectionYdummy->GetMeanError());
      h_JER_pT[iInp][irad]->SetBinContent(i,projectionYdummy->GetStdDev()); //GetMeanError()
      h_JER_pT[iInp][irad]->SetBinError(i,projectionYdummy->GetRMSError()); //GetMeanError()
    }
  }
}

// 2D PLOT
  for(Int_t irad=0; irad<3;irad++){
    TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
    DrawGammaCanvasSettings( cReso, 0.11, 0.01, 0.01, 0.105);
    // cReso->SetLogz();

    TH2F* histoJESDummy   = new TH2F("histoJESDummy","histoJESDummy",1000,2.5, 29.9,1000,-0.5, 0.5);
    SetStyleHistoTH2ForGraphs(histoJESDummy, "#it{p}_{T}^{true} (GeV/#it{c})","(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.99);
    histoJESDummy->GetXaxis()->SetNoExponent();
    histoJESDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histoJESDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoJESDummy->Draw();
    TLegend* legendJES  = GetAndSetLegend2(0.13, 0.83-(4*textSizeLabelsRel), 0.6, 0.83-(1*textSizeLabelsRel),1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    DrawGammaLines(0, 29, 0., 0., 2, kGray+2, 7);
    for(Int_t iInp=0; iInp<3;iInp++){
      DrawGammaSetMarker(h_JES_pT[iInp][irad], markerInput[iInp], 3, colorInput[iInp], colorInput[iInp]);
      h_JES_pT[iInp][irad]->Draw("same,p");
      legendJES->AddEntry(h_JES_pT[iInp][irad],Form("%s",namePlotInput[iInp].Data()),"p");
    }
    legendJES->Draw();
    drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot[irad],etaMin[irad],etaMax[irad]),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cReso->Print(Form("%s/JES_jet_pT_%s.%s", outputDir.Data(),jetradius[irad].Data(), suffix.Data()));


    TH2F* histoJERDummy   = new TH2F("histoJERDummy","histoJERDummy",1000,2.5, 29.9,1000,-0.01, 0.39);
    SetStyleHistoTH2ForGraphs(histoJERDummy, "#it{p}_{T}^{true} (GeV/#it{c})","#sigma(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.99);
    histoJERDummy->GetXaxis()->SetNoExponent();
    histoJERDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histoJERDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoJERDummy->Draw();
    for(Int_t iInp=0; iInp<3;iInp++){
      DrawGammaSetMarker(h_JER_pT[iInp][irad], markerInput[iInp], 3, colorInput[iInp], colorInput[iInp]);
      h_JER_pT[iInp][irad]->Draw("same,p");
      // drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      // drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }
    legendJES->Draw();
       DrawGammaLines(2.5,29.9,0.1,0.1,2,kGray+1,7);
   drawLatexAdd(collisionSystem.Data(),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("Anti-#it{k}_{T}, HCAL+ECAL jets"),0.15,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#it{R} < %1.1f, %1.1f < #eta_{jet} < %1.1f",jetradiusplot[irad],etaMin[irad],etaMax[irad]),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cReso->Print(Form("%s/JER_jet_pT_%s.%s", outputDir.Data(), jetradius[irad].Data(),suffix.Data()));

  }
  for(Int_t iInp=0; iInp<3;iInp++){

    TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
    DrawGammaCanvasSettings( cReso, 0.1, 0.01, 0.01, 0.105);

    TH2F * histResoDummy   = new TH2F("histResoDummy","histResoDummy",1000,2.5, 29.9,1000,-0.5, 0.5);
    histResoDummy->GetXaxis()->SetNoExponent();
    histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
    histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    histResoDummy->GetZaxis()->SetRangeUser(0,500);
    for(Int_t irad=0; irad<3;irad++){
      // histResoDummy->Draw();
      SetStyleHistoTH2ForGraphs(h_jetscale_pT[iInp][irad], "#it{p}_{T}^{true} (GeV/#it{c})","(#it{p}_{T}^{rec} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);

      h_jetscale_pT[iInp][irad]->Draw("col");
      DrawGammaSetMarker(h_JES_pT[iInp][irad], 20, 2, kMagenta+1, kMagenta+1);
      h_JES_pT[iInp][irad]->Draw("same");
      drawLatexAdd(collisionSystem.Data(),0.90,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s",namePlotInput[iInp].Data()),0.90,0.70,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/2D_jet_pT_%s_%s.%s", outputDir.Data(), inputFolderNames[iInp].Data(),jetradius[irad].Data(), suffix.Data()));
    }
    for(Int_t irad=0; irad<3;irad++){
      // histResoDummy->Draw();
      SetStyleHistoTH2ForGraphs(h_jet_dEtadPhi[iInp][irad], "#eta^{rec}-#eta^{true}","#phi^{rec}-#phi^{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);

      h_jet_dEtadPhi[iInp][irad]->Draw("col");
      drawLatexAdd(collisionSystem.Data(),0.90,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("Anti-#it{k}_{T}, #it{R} = %1.1f",jetradiusplot[irad]),0.90,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%1.1f < #eta_{jet} < %1.1f",etaMin[irad],etaMax[irad]),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s",namePlotInput[iInp].Data()),0.90,0.70,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      DrawGammaLines(h_jet_dEtadPhi[iInp][irad]->GetMean(1), h_jet_dEtadPhi[iInp][irad]->GetMean(1), -0.5, 0.5, 2, kGray+2, 7);
      DrawGammaLines( -0.5, 0.5,h_jet_dEtadPhi[iInp][irad]->GetMean(2), h_jet_dEtadPhi[iInp][irad]->GetMean(2), 2, kGray+2, 7);
     cReso->Print(Form("%s/dEtadPhi_%s_%s.%s", outputDir.Data(), inputFolderNames[iInp].Data(),jetradius[irad].Data(), suffix.Data()));
    }
    // TH2F * histResoDummyE   = new TH2F("histResoDummyE","histResoDummyE",1000,0, 119,1000,-0.5, 0.5);
    // SetStyleHistoTH2ForGraphs(histResoDummyE, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    // histResoDummyE->GetXaxis()->SetNoExponent();
    // histResoDummyE->GetYaxis()->SetNdivisions(505,kTRUE);
    // histResoDummyE->GetXaxis()->SetMoreLogLabels(kTRUE);
    // histResoDummyE->Draw();

    // h_jetscale_E[irad]->Draw("col,same");
    // cReso->Print(Form("%s/2D_jet_E_%s.%s", outputDir.Data(), inputFolderNames[iInp].Data(), suffix.Data()));

  }


}