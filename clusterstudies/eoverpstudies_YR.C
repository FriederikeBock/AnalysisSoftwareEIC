#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 

float getEPcutvalue(int icalo, float partE, bool useFWHM = false){
  float EP_cut_strength = 1.6;
  // float EP_cut_strength = 4.0;
  float EP_cut = 1;
  if(useFWHM){
    if(icalo==0) EP_cut = 1. - EP_cut_strength * std::sqrt( 2.1 * 2.1 + 1.4 * 1.4 / partE ) / 100.;
    if(icalo==1) EP_cut = 1. - EP_cut_strength * std::sqrt( 1.6 * 1.6 + 1.3 * 1.3 / partE ) / 100.;
    if(icalo==2) EP_cut = 1. - EP_cut_strength * std::sqrt( 0.3 * 0.3 + 7.8 * 7.8 / partE ) / 100.;
  } else {
    if(icalo==0) EP_cut = 1. - EP_cut_strength * std::sqrt( 0.8 * 0.8 + 2.0 * 2.0 / partE ) / 100.;
    if(icalo==1) EP_cut = 1. - EP_cut_strength * std::sqrt( 1.3 * 1.3 + 0.8 * 0.8 / partE ) / 100.;
    if(icalo==2) EP_cut = 1. - EP_cut_strength * std::sqrt( 0.0 * 0.0 + 7.3 * 7.3 / partE ) / 100.;
  }
  return EP_cut;
}
void eoverpstudies_YR(
    TString inputFileNamePion     = "central_allculsters.root",
    TString inputFileNameElectron     = "central_allculsters.root",
    TString suffix            = "pdf",
    Bool_t doFitting          = kTRUE,
    TString collisionsSys     = "PythiaMB"
){  
  TString clusterizerName = "MA clusters";
  const int maxcalo = 3;
  TString caloName[maxcalo+1] = {"EEMC", "BECAL", "FEMC", "Hybrid"};
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
  TString perfLabel                 = "ECCE G4 simulation";
  

  Color_t colorPID_light[nPID]          = {kGray+1, kRed-6, kGreen-6, kCyan-6, kBlue+1, kOrange};
  Color_t colorPID_light2[nPID]         = {kGray+2, kRed-2, kGreen-2, kCyan-6, kBlue+1, kOrange};

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
  
  TString outputDir                 = Form("plotsEoverP");
  if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFilePion                      = new TFile(inputFileNamePion.Data());
  if(inputFilePion->IsZombie()){
    cout << "inputFilePion not found!" << endl;
  }
  TFile* inputFileElec                      = new TFile(inputFileNameElectron.Data());
  if(inputFileElec->IsZombie()){
    cout << "inputFileElec not found!" << endl;
  }

  //************************** Read data **************************************************
  const int nEne = 15;
  const static Double_t partE[]   = {   0.7, 1.0, 1.5, 2, 3, 4, 6, 10,15,20,25,  30,35,  40, 50};
  // double energymax = 40.;
  // const int nEne = 26;
  // const static Double_t partE[]   = {   0.7, 0.9, 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 4.5,  5,
  //                                       5.5, 6, 7, 8, 10,
  //                                      12.5,15,17.5,20,22.5,25,  30,35,  40, 50};
  // double energymax = 40.;
  int padebin = 0;
  int padebinelec = 0;
  const int maxmatch = 2;
  TH1F*  h_EoverP_trackSpec_pions[maxcalo] = {NULL}; // [calorimeter_enum][algorithm_enum]
  TH1F*  h_EoverP_trackSpec_electrons[maxcalo] = {NULL}; // [calorimeter_enum][algorithm_enum]
  TH2F*  h_EoverP_pions_EoverP_E[maxmatch][maxcalo] = {{NULL}}; // [calorimeter_enum][algorithm_enum]
  TH2F*  h_EoverP_electrons_EoverP_E[maxmatch][maxcalo] = {{NULL}}; // [calorimeter_enum][algorithm_enum]
  TH3F*  h_EoverP_pions_EoverP_E_Eta[maxmatch][maxcalo][Nusefilltrue] = {{{NULL}}}; // [calorimeter_enum][algorithm_enum]
  TH3F*  h_EoverP_electrons_EoverP_E_Eta[maxmatch][maxcalo][Nusefilltrue] = {{{NULL}}}; // [calorimeter_enum][algorithm_enum]

  TH1F* histEoP_proj_pions[Nusefilltrue][maxmatch][maxcalo][nEta+1][50]     = {{{{{NULL}}}}};
  TH1F* histEoP_proj_electrons[Nusefilltrue][maxmatch][maxcalo][nEta+1][50] = {{{{{NULL}}}}};
  TH1F* histIntCounts_pions_EoP[Nusefilltrue][maxmatch][maxcalo][nEta+1]    = {{{{NULL}}}};
  TH1F* histIntCounts_pions_EoP_fwhm[Nusefilltrue][maxmatch][maxcalo][nEta+1]    = {{{{NULL}}}};
  TH1F* histIntCounts_pions_EoP_95[Nusefilltrue][maxmatch][maxcalo][nEta+1]    = {{{{NULL}}}};
  TH1F* histIntCounts_electrons_EoP[Nusefilltrue][maxmatch][maxcalo][nEta+1]= {{{{NULL}}}};
  TH1F* histIntCounts_electrons_EoP_fwhm[Nusefilltrue][maxmatch][maxcalo][nEta+1]= {{{{NULL}}}};
  TH1F* histIntCounts_pions_all[Nusefilltrue][maxmatch][maxcalo][nEta+1]    = {{{{NULL}}}};
  TH1F* histIntCounts_electrons_all[Nusefilltrue][maxmatch][maxcalo][nEta+1]= {{{{NULL}}}};
  TH1F* histPionRejection[Nusefilltrue][maxmatch][maxcalo][nEta+1]          = {{{{NULL}}}};
  TH1F* histPionRejection_fwhm[Nusefilltrue][maxmatch][maxcalo][nEta+1]          = {{{{NULL}}}};
  TH1F* histPionRejection_95[Nusefilltrue][maxmatch][maxcalo][nEta+1]          = {{{{NULL}}}};
  TH1F* histElectronEfficiency[Nusefilltrue][maxmatch][maxcalo][nEta+1]          = {{{{NULL}}}};
  TH1F* histElectronEfficiency_fwhm[Nusefilltrue][maxmatch][maxcalo][nEta+1]          = {{{{NULL}}}};

  double EP_cut[Nusefilltrue][maxmatch][maxcalo][nEta+1][50] = {0};
  double EP_cut_95[Nusefilltrue][maxmatch][maxcalo][nEta+1][50] = {0};
  double EP_cut_FWHM[Nusefilltrue][maxmatch][maxcalo][nEta+1][50] = {0};
  TF1* EP_fit[Nusefilltrue][maxmatch][maxcalo][nEta+1][50] = {NULL};

  TH2F * histProjDummy                           = new TH2F("histProjDummy","histProjDummy",1000,0, 2,1000,0.8, 1e6);
  SetStyleHistoTH2ForGraphs(histProjDummy, "#it{p}_{track}", "counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  // TString matchbasis = "trackbased";
  TString matchbasis[maxmatch] = {"trackbased", "projectionbased"};
  int imatch=1;
  for(int ifilltrue=0;ifilltrue<Nusefilltrue;ifilltrue++){
  // for(int imatch=0;imatch<maxmatch;imatch++){
    for(int icalo=0;icalo<maxcalo;icalo++){

      gSystem->Exec(Form("mkdir -p %s/%s/%s",outputDir.Data(), caloName[icalo].Data(), matchbasis[0].Data()));
      gSystem->Exec(Form("mkdir -p %s/%s/%s",outputDir.Data(), caloName[icalo].Data(), matchbasis[1].Data()));
      cout << "processing " << caloName[icalo].Data() << endl;
      TString NameEP_hist_EpE = Form("%s/h_EoverP_pions_EoverP_E_%s_%s",caloName[icalo].Data(),matchbasis[imatch].Data(), caloName[icalo].Data() );
      // TString NameE_Minv = Form("%s/h_Etrue_Minv_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
      h_EoverP_pions_EoverP_E[imatch][icalo]    = (TH2F*)inputFilePion->Get(NameEP_hist_EpE.Data());
      if(!h_EoverP_pions_EoverP_E[imatch][icalo]){
        cout << "\t" <<NameEP_hist_EpE.Data() << " not found" << endl;
        continue;
      } else {
        caloactive[icalo]=true;
      };
      h_EoverP_electrons_EoverP_E[imatch][icalo]    = (TH2F*)inputFileElec->Get(Form("%s/h_EoverP_electrons_EoverP_E_%s_%s",caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data() ));
      h_EoverP_trackSpec_pions[icalo]    = (TH1F*)inputFilePion->Get(Form("%s/h_EoverP_trackSpec_pions_%s",caloName[icalo].Data(), caloName[icalo].Data() ));
      h_EoverP_trackSpec_electrons[icalo]    = (TH1F*)inputFileElec->Get(Form("%s/h_EoverP_trackSpec_electrons_%s",caloName[icalo].Data(), caloName[icalo].Data() ));

      h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]    = (TH3F*)inputFilePion->Get(Form("%s/h_EoverP_pions_EoverP_E_Eta_%s_%s_%s",caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data() ));
      h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]    = (TH3F*)inputFileElec->Get(Form("%s/h_EoverP_electrons_EoverP_E_Eta_%s_%s_%s",caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data() ));




      if (h_EoverP_pions_EoverP_E[imatch][icalo]->GetEntries() > 100 || h_EoverP_electrons_EoverP_E[imatch][icalo]->GetEntries()>100){
        float latextopedge = 0.90;
        float latexsideedge = 0.92;
        float  textSizeLabelsRel      = 55./1200;
        TCanvas* Minvpi02D = new TCanvas(Form("Minvpi02D%s",caloName[icalo].Data()),"",0,0,1500,800);
        DrawGammaCanvasSettings( Minvpi02D, 0.095, 0.01, 0.01, 0.105);
        Minvpi02D->Divide(2,1);
        Minvpi02D->cd(1);
        h_EoverP_pions_EoverP_E[imatch][icalo]->GetXaxis()->SetRangeUser(0.,2);
        // h_EoverP_pions_EoverP_E[imatch][icalo]->GetYaxis()->SetRangeUser(0.,0.2);
        SetStyleHistoTH2ForGraphs(h_EoverP_pions_EoverP_E[imatch][icalo],"E^{MC}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
        h_EoverP_pions_EoverP_E[imatch][icalo]->Draw("colz");
        drawLatexAdd("true #pi^{-} tracks",0.17,latextopedge,textSizeLabelsRel,false,false,false);      
        Minvpi02D->cd(2);
        h_EoverP_electrons_EoverP_E[imatch][icalo]->GetXaxis()->SetRangeUser(0.,2);
        // h_EoverP_electrons_EoverP_E[imatch][icalo]->GetYaxis()->SetRangeUser(0.,0.2);
        SetStyleHistoTH2ForGraphs(h_EoverP_electrons_EoverP_E[imatch][icalo],"E^{reco}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
        h_EoverP_electrons_EoverP_E[imatch][icalo]->Draw("colz");
        
        // TLatex Tl1;
        // Tl1.SetTextSize(0.035);
        // Tl1.DrawLatex(0.55*energymax,0.04,perfLabel);
        // Tl1.DrawLatex(0.55*energymax,0.03,Form("Calorimeter: %s ", caloName[icalo].Data()));
        // Tl1.DrawLatex(0.55*energymax,0.02,Form("Clusterization type: %s ", clusterizerName.Data()));
        // Tl1.DrawLatex(0.55*energymax,0.01,"#pi ^{0} #rightarrow #gamma #gamma");
        // Tl1.Draw();

        drawLatexAdd(perfLabel.Data(),latexsideedge,latextopedge,textSizeLabelsRel,false,false,true);
        drawLatexAdd(collisionSystem.Data(),latexsideedge,latextopedge-0.05,textSizeLabelsRel,false,false,true);
        drawLatexAdd(Form("Calorimeter: %s", caloName[icalo].Data()),latexsideedge,latextopedge-0.1,textSizeLabelsRel,false,false,true);      
        drawLatexAdd("true e^{-} tracks",0.17,latextopedge,textSizeLabelsRel,false,false,false);      
        Minvpi02D->Print(Form("%s/%s/%s/EoverP_vs_E_%s_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(), str_EoverP_FillTrue[ifilltrue].Data(),suffix.Data()));
      }

      float minEtaProj = 0;
      float maxEtaProj = 0;

      for(int ieta=0;ieta<nEta+1;ieta++){

        if(ieta==nEta){
            minEtaProj = nominalEtaRegion[icalo][0];
            maxEtaProj = nominalEtaRegion[icalo][1];
        } else {
          minEtaProj = partEtaCalo[ieta];
          maxEtaProj = partEtaCalo[ieta+1];
        }
        cout << caloName[icalo] << endl;
        cout << minEtaProj << " < eta < " << maxEtaProj << endl;
        padebin = 0;
        padebinelec = 0;
        histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_pions_EoP_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_pions_EoP_fwhm[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_pions_EoP_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_pions_EoP_95[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_pions_EoP_95_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_electrons_EoP_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_electrons_EoP_fwhm[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_electrons_EoP_fwhm_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_pions_all_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_electrons_all_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histPionRejection[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histPionRejection_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histPionRejection_fwhm[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histPionRejection_fwhm_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histPionRejection_95[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histPionRejection_95_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histElectronEfficiency[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histElectronEfficiency_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histElectronEfficiency_fwhm[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histElectronEfficiency_fwhm_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);

        for(int ebin=0;ebin<nEne-1;ebin++){

          h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->SetRange(  h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin]),h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin+1]));
          h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->SetRange(  h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(minEtaProj),h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(maxEtaProj));
          h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->SetRange( h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin]),h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin+1]));
          h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->SetRange( h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(minEtaProj),h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(maxEtaProj));

          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("%1.1f_%1.1f_%d_%s_x",minEtaProj,partE[ebin],imatch,caloName[icalo].Data()));
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("%1.1f_%1.1f_%d_%s_x",minEtaProj,partE[ebin],imatch,caloName[icalo].Data()));
          // histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("histEoP_proj_pions_%s_%f_%f",caloName[icalo].Data(),minEtaProj, partE[ebin]));
          // histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("histEoP_proj_pions_%s_%f_%f",caloName[icalo].Data(),minEtaProj, partE[ebin]));
          float EP_cut_tmp = getEPcutvalue(icalo,(partE[ebin]+partE[ebin+1])/2);
          float EP_cut_tmp_fwhm = getEPcutvalue(icalo,(partE[ebin]+partE[ebin+1])/2,true);
          double EP_shift =0;
          if(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetEntries()>10){
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]     = new TF1(Form("fit_%1.1f_%1.1f_%d_%s_x",minEtaProj,partE[ebin],imatch,caloName[icalo].Data()), "gaus", 0.8, 1.05);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParLimits(1, 0.8, 1.0);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParameter(1, 0.95);
            double fitsigmatemp = 0.02;
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParLimits(2, 0.6*fitsigmatemp,3*fitsigmatemp);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParameter(2, fitsigmatemp);
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->Fit(EP_fit[ifilltrue][imatch][icalo][ieta][padebin],"L0RMEQ","",0.8, 1.05);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParLimits(1, EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1)*0.9, EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1)*1.1);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParameter(1, EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1));
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParLimits(2, EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(2)*0.5,EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(2)*1.5);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetParameter(2, EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(2));
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->Fit(EP_fit[ifilltrue][imatch][icalo][ieta][padebin],"L0RMEQ","",EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1)-EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(2)*2, EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1)+EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(2)*3);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->SetRange(EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1)-EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(2)*2, EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1)+EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(2)*3);
            // EP_shift     = 1 - EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1);
            double EP_shift_bin = FindLargestBin1DHist(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin], 0.6, 1.1 );
            // cout << EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParameter(1) << " vs " << EP_shift_bin << endl;
            EP_shift     = 1 - EP_shift_bin;
            // meanErr  = EP_fit[ifilltrue][imatch][icalo][ieta][padebin]->GetParError(1);
            // EP_shift = 1-histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetMaximum();
            EP_cut[ifilltrue][imatch][icalo][ieta][padebin] = EP_cut_tmp-EP_shift;
            EP_cut_FWHM[ifilltrue][imatch][icalo][ieta][padebin] = EP_cut_tmp_fwhm-EP_shift;
          }else {
            // if(padebin==0)
              EP_cut[ifilltrue][imatch][icalo][ieta][padebin] = EP_cut_tmp;
              EP_cut_FWHM[ifilltrue][imatch][icalo][ieta][padebin] = EP_cut_tmp_fwhm;
              EP_fit[ifilltrue][imatch][icalo][ieta][padebin]     = new TF1(Form("fit_%1.1f_%1.1f_%d_%s_x",minEtaProj,partE[ebin],imatch,caloName[icalo].Data()), "gaus", 0.8, 1.05);
            // else
            //   EP_cut[ifilltrue][imatch][icalo][ieta][padebin] = EP_cut[ifilltrue][imatch][icalo][ieta][padebin-1];
          }
          // cout << "EPcutvalue: " << EP_cut[ifilltrue][imatch][icalo][ieta][padebin] << endl;

          Double_t intpions  = histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(EP_cut[ifilltrue][imatch][icalo][ieta][padebin]),
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          if(intpions==0)intpions=1;
          histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpions);

          Double_t intpionsall  = histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(0.0),
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpionsall);
          // padebin++;

          Double_t intelectrons  = histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(EP_cut[ifilltrue][imatch][icalo][ieta][padebin]),
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          if(intelectrons==0)intelectrons=1;
          histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectrons);

          Double_t intelectronsall  = histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(0.0),
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          if(intelectronsall==0)intelectronsall=1;
          histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectronsall);

          Double_t intelectrons_fwhm  = histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(EP_cut_FWHM[ifilltrue][imatch][icalo][ieta][padebin]),
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          if(intelectrons_fwhm==0)intelectrons_fwhm=1;
          histIntCounts_electrons_EoP_fwhm[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_electrons_EoP_fwhm[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectrons_fwhm);

          Double_t intelectrons_95 = intelectrons;
          EP_cut_95[ifilltrue][imatch][icalo][ieta][padebin] = EP_cut[ifilltrue][imatch][icalo][ieta][padebin];
          for(int iwindow=0;iwindow<100;iwindow++){
            if(abs(intelectrons_95/intelectronsall)>0.95) break;
            else {
              intelectrons_95 = histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
              histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(EP_cut[ifilltrue][imatch][icalo][ieta][padebin]-0.005*iwindow),
              histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
              EP_cut_95[ifilltrue][imatch][icalo][ieta][padebin] = EP_cut[ifilltrue][imatch][icalo][ieta][padebin]-0.005*iwindow;
            }
          }


          Double_t intpions_95  = histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(EP_cut_95[ifilltrue][imatch][icalo][ieta][padebin]),
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          if(intpions_95==0)intpions_95=1;
          histIntCounts_pions_EoP_95[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_EoP_95[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpions_95);


          Double_t intpions_fwhm  = histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->Integral(
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(EP_cut_FWHM[ifilltrue][imatch][icalo][ieta][padebin]),
            histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          if(intpions_fwhm==0)intpions_fwhm=1;
          histIntCounts_pions_EoP_fwhm[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_EoP_fwhm[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpions_fwhm);


          // cout << (partE[ebin]+partE[ebin+1])/2 << "\tEP-default: " << EP_cut[ifilltrue][imatch][icalo][ieta][padebin] << "\tEP-95: " << EP_cut_95[ifilltrue][imatch][icalo][ieta][padebin] << "\tsurvival: " << intelectrons_95/intelectronsall << endl;
          double rejection = histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) / histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2));
          if(isnan(rejection)) rejection = -1;
          if(isinf(rejection)) rejection = -1;
          // histPionRejection[ifilltrue][imatch][icalo][ieta]->SetBinContent(histPionRejection[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),rejection);
          // histElectronEfficiency[ifilltrue][imatch][icalo][ieta]->SetBinContent(histElectronEfficiency[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectrons/intelectronsall);
          // histElectronEfficiency[ifilltrue][imatch][icalo][ieta]->SetBinError(histElectronEfficiency[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),0.03*intelectrons/intelectronsall);
          padebin++;
        }
      histPionRejection[ifilltrue][imatch][icalo][ieta]->Divide(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta],histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta], 1, 1,"B" );
      histPionRejection_fwhm[ifilltrue][imatch][icalo][ieta]->Divide(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta],histIntCounts_pions_EoP_fwhm[ifilltrue][imatch][icalo][ieta], 1, 1,"B" );
      histPionRejection_95[ifilltrue][imatch][icalo][ieta]->Divide(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta],histIntCounts_pions_EoP_95[ifilltrue][imatch][icalo][ieta], 1, 1,"B" );
      histElectronEfficiency[ifilltrue][imatch][icalo][ieta]->Divide(histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta],histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta], 1, 1,"B" );
      histElectronEfficiency_fwhm[ifilltrue][imatch][icalo][ieta]->Divide(histIntCounts_electrons_EoP_fwhm[ifilltrue][imatch][icalo][ieta],histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta], 1, 1,"B" );
      TCanvas* mass2D = new TCanvas("mass2D","",0,0,1000,800);
      gStyle->SetTitleX(0.5);
      gStyle->SetTitleAlign(23);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleFontSize(0.07); 
      DrawGammaCanvasSettings( mass2D, 0.095, 0.01, 0.01, 0.10);
      mass2D->Divide(5,6);
      int padnumber = 0;
      TLegend* legend1 = new TLegend(0., 0.01, 0.4, 0.25);
      for(int z=0; z< padebin +1 ; z++ ){
        if(z==4){
          mass2D->cd(z+1);
          TLatex Tl;
          Tl.SetTextSize(0.1);
          Tl.DrawLatex(0.03,0.9,perfLabel);
          Tl.DrawLatex(0.03,0.75,Form("Calorimeter: %s ", caloName[icalo].Data()));
          Tl.DrawLatex(0.03,0.6,Form("Clusterization type: %s ", clusterizerName.Data()));
          // Tl.DrawLatex(0.03,0.45,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.45,Form("%1.2f < #eta < %1.1f", minEtaProj, maxEtaProj));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          gPad->SetLogy();
          if(!histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          histProjDummy->GetYaxis()->SetRangeUser(0.8,histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum()*10);
          histProjDummy->DrawCopy();
          // SetStyleHistoTH1ForGraphs( histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
          DrawGammaSetMarker(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber],34,1, kBlack, kBlack);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->Draw("same,ep");

        drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        padnumber++;
        }
      }
      if(minEtaProj >= etarangemin[icalo] && maxEtaProj <= etarangemax[icalo])
        mass2D->Print(Form("%s/%s/%s/EoP_bin_pions%s_%s_%2.1f_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(), minEtaProj,str_EoverP_FillTrue[ifilltrue].Data(), suffix.Data()));

      mass2D = new TCanvas("mass2D","",0,0,1000,800);
      gStyle->SetTitleX(0.5);
      gStyle->SetTitleAlign(23);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleFontSize(0.07); 
      DrawGammaCanvasSettings( mass2D, 0.095, 0.01, 0.01, 0.10);
      mass2D->Divide(5,6);
      padnumber = 0;
      legend1 = new TLegend(0., 0.01, 0.4, 0.25);
      for(int z=0; z< padebin +1 ; z++ ){
        if(z==4){
          mass2D->cd(z+1);
          TLatex Tl;
          Tl.SetTextSize(0.1);
          Tl.DrawLatex(0.03,0.9,perfLabel);
          Tl.DrawLatex(0.03,0.75,Form("Calorimeter: %s ", caloName[icalo].Data()));
          Tl.DrawLatex(0.03,0.6,Form("Clusterization type: %s ", clusterizerName.Data()));
          // Tl.DrawLatex(0.03,0.45,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.45,Form("%1.2f < #eta < %1.1f", minEtaProj, maxEtaProj));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          gPad->SetLogy();
          if(!histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          histProjDummy->GetYaxis()->SetRangeUser(0.8,histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum()*10);
          histProjDummy->DrawCopy();
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
          DrawGammaSetMarker(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber],34,1, kBlack, kBlack);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->Draw("same,ep");
          float EP_cut_tmp = EP_cut[ifilltrue][imatch][icalo][ieta][padnumber];
            DrawGammaLines(EP_cut_tmp, EP_cut_tmp, 0, 0.5*histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            // DrawGammaLines(eopcutvaluemax[icalo], eopcutvaluemax[icalo], 0, 0.5*histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
          //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
          //   if(z==0){
          //     legend1->AddEntry(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber], "Simulation", "p");
          //     legend1->AddEntry(fithistEtrueMinvFINAL[icalo][padnumber], "Fit", "l");
          //   }
        drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        padnumber++;
        }
      }

      if(minEtaProj >= etarangemin[icalo] && maxEtaProj <= etarangemax[icalo])
        mass2D->Print(Form("%s/%s/%s/EoP_bin_electrons%s_%s_%2.1f_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(), minEtaProj,str_EoverP_FillTrue[ifilltrue].Data(), suffix.Data()));
      if(ieta==nEta)
        mass2D->Print(Form("%s/%s/%s/EoP_bin_electrons%s_%s_full_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(),str_EoverP_FillTrue[ifilltrue].Data(), suffix.Data()));
      mass2D = new TCanvas("mass2D","",0,0,1000,800);
      gStyle->SetTitleX(0.5);
      gStyle->SetTitleAlign(23);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleFontSize(0.07); 
      DrawGammaCanvasSettings( mass2D, 0.095, 0.01, 0.01, 0.10);
      mass2D->Divide(5,6);
      padnumber = 0;
      legend1 = new TLegend(0., 0.01, 0.4, 0.25);
      for(int z=0; z< padebin +1 ; z++ ){
        if(z==4){
          mass2D->cd(z+1);
          TLatex Tl;
          Tl.SetTextSize(0.1);
          Tl.DrawLatex(0.03,0.9,perfLabel);
          Tl.DrawLatex(0.03,0.75,Form("Calorimeter: %s ", caloName[icalo].Data()));
          Tl.DrawLatex(0.03,0.6,Form("Clusterization type: %s ", clusterizerName.Data()));
          // Tl.DrawLatex(0.03,0.45,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.45,Form("%1.2f < #eta < %1.1f", minEtaProj, maxEtaProj));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          gPad->SetLogy();
          if(!histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          histProjDummy->GetYaxis()->SetRangeUser(0.8,histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum()*10);
          histProjDummy->DrawCopy();
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber], "#it{p}_{track}", "counts", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,0.6,1.1);
          DrawGammaSetMarker(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber],34,1, kRed+2, kRed+2);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->Draw("same,ep");
          float EP_cut_tmp = EP_cut[ifilltrue][imatch][icalo][ieta][padnumber];
            DrawGammaLines(EP_cut_tmp, EP_cut_tmp, 0, 0.5*histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            // DrawGammaLines(eopcutvalue[icalo], eopcutvalue[icalo], 0, 0.5*histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            // DrawGammaLines(eopcutvaluemax[icalo], eopcutvaluemax[icalo], 0, 0.5*histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
          //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
            if(z==0){
              legend1->AddEntry(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber], "Pion E/p", "p");
            }
        drawLatexAdd(Form("%1.1f < #it{p}_{track} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        padnumber++;
        }
      }

      if(minEtaProj >= etarangemin[icalo] && maxEtaProj <= etarangemax[icalo])
      mass2D->Print(Form("%s/%s/%s/EoP_bin_pions_%s_%s_%2.1f_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(), minEtaProj,str_EoverP_FillTrue[ifilltrue].Data(), suffix.Data()));
      if(ieta==nEta)
      mass2D->Print(Form("%s/%s/%s/EoP_bin_pions_%s_%s_full_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(),str_EoverP_FillTrue[ifilltrue].Data(), suffix.Data()));

      mass2D = new TCanvas("mass2D","",0,0,2000,1600);
      gStyle->SetTitleX(0.5);
      gStyle->SetTitleAlign(23);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleFontSize(0.07); 
      DrawGammaCanvasSettings( mass2D, 0.095, 0.01, 0.01, 0.10);
      mass2D->Divide(5,6);
      padnumber = 0;
      legend1 = new TLegend(0., 0.01, 0.4, 0.25);
      for(int z=0; z< padebin +1 ; z++ ){
        if(z==4){
          mass2D->cd(z+1);
          TLatex Tl;
          Tl.SetTextSize(0.1);
          Tl.DrawLatex(0.03,0.9,perfLabel);
          Tl.DrawLatex(0.03,0.75,Form("Calorimeter: %s ", caloName[icalo].Data()));
          Tl.DrawLatex(0.03,0.6,Form("Clusterization type: %s ", clusterizerName.Data()));
          // Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.45,Form("%1.2f < #eta < %1.1f", minEtaProj, maxEtaProj));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          gPad->SetLogy();
          if(!histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          histProjDummy->GetYaxis()->SetRangeUser(0.8,histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum()*10);
          histProjDummy->DrawCopy();
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          // SetStyleHistoTH1ForGraphs( histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber], "#it{p}_{track}", "counts", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,0.6,1.1);
          DrawGammaSetMarker(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber],34,1, kBlack, kBlack);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->Draw("same,ep");
          if(EP_fit[ifilltrue][imatch][icalo][ieta][padnumber]){
          DrawGammaSetMarkerTF1(EP_fit[ifilltrue][imatch][icalo][ieta][padnumber], 1, 1, kBlue+2);
          // EP_fit[ifilltrue][imatch][icalo][ieta][padnumber]->SetRange(minEPart[i], maxEPart[i]);
          EP_fit[ifilltrue][imatch][icalo][ieta][padnumber]->Draw("same");
          }
          DrawGammaSetMarker(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber],28,1, kRed+2, kRed+2);
          // histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          // histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          // histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber]->Draw("ep,same");
          // DrawGammaSetMarker(histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
          // histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
          // histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber]->Draw("same,hist");
          // if(fithistEtrueMinvFINAL[icalo][padnumber]){
          //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
          //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
          //   fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
          
          float EP_cut_tmp = EP_cut[ifilltrue][imatch][icalo][ieta][padnumber];
            DrawGammaLines(EP_cut_tmp, EP_cut_tmp, 0, 0.5*histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            // DrawGammaLines(eopcutvalue[icalo], eopcutvalue[icalo], 0, 0.5*histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            // DrawGammaLines(eopcutvaluemax[icalo], eopcutvaluemax[icalo], 0, 0.5*histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
          //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
            if(z==0){
              legend1->AddEntry(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padnumber], "Electron E/p", "p");
              legend1->AddEntry(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padnumber], "Pion E/p", "p");
            }
        drawLatexAdd(Form("%1.1f < #it{p}_{track} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        padnumber++;
        }
      }

      if(minEtaProj >= etarangemin[icalo] && maxEtaProj <= etarangemax[icalo])
        mass2D->Print(Form("%s/%s/%s/EoP_bin_electrons_and_pions_%s_%s_%2.1f_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(), minEtaProj,str_EoverP_FillTrue[ifilltrue].Data(), suffix.Data()));
      if(ieta==nEta)
        mass2D->Print(Form("%s/%s/%s/EoP_bin_electrons_and_pions_%s_%s_full_%s.%s", outputDir.Data(), caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data(), clusterizerName.Data(),str_EoverP_FillTrue[ifilltrue].Data(), suffix.Data()));
    }
  }
  }

cout << "meow" << endl;
//   // **********************************************************************************************************************
//   // ******************************** Acceptance * Efficiency for pi0 single measurement 8TeV **************************
//   // **********************************************************************************************************************
  textSizeLabelsPixel             = 55;
  textSizeLabelsRel      = 55./1200;
  cout << textSizeLabelsRel << endl;

  TCanvas* canvasMergingFraction       = new TCanvas("canvasMergingFraction", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasMergingFraction,  0.1, 0.01, 0.015, 0.095);
  canvasMergingFraction->SetLogy(1);
  canvasMergingFraction->SetLogx(1);

  TH2F * histo2DAccEff;
  histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.51, 49, 1000, 1, 99999 );
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{track} (GeV/#it{c})", "#pi^{#pm}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
  histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
  histo2DAccEff->GetXaxis()->SetNoExponent();
  histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
  TGraphAsymmErrors* graphPionRejection[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
  TGraphAsymmErrors* graphPionRejection_fwhm[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
  TGraphAsymmErrors* graphPionRejection_95[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
  for(int ifilltrue=0;ifilltrue<Nusefilltrue;ifilltrue++){
    SetStyleHistoTH2ForGraphs( histo2DAccEff, ifilltrue==0 ? "#it{p}_{track, rec} (GeV/#it{c})" :  "#it{p}_{track, true} (GeV/#it{c})", "#pi^{#pm}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
    for (Int_t icalo = 0; icalo < maxcalo; icalo++){
      cout << caloName[icalo].Data() << endl;
      histo2DAccEff->DrawCopy();
      for (Int_t ieta = 0; ieta < nEta+1; ieta++){
        if(ieta==nEta){
          cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
          if(histPionRejection[ifilltrue][imatch][icalo][ieta]){
            graphPionRejection[ifilltrue][icalo][ieta] = new TGraphAsymmErrors(histPionRejection[ifilltrue][imatch][icalo][ieta]);
            while(graphPionRejection[ifilltrue][icalo][ieta]->GetN()>0 && graphPionRejection[ifilltrue][icalo][ieta]->GetY()[graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1] == -1) graphPionRejection[ifilltrue][icalo][ieta]->RemovePoint(graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1);
            while(graphPionRejection[ifilltrue][icalo][ieta]->GetN()>0 && graphPionRejection[ifilltrue][icalo][ieta]->GetX()[graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1] > energymax[icalo]) graphPionRejection[ifilltrue][icalo][ieta]->RemovePoint(graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1);
            DrawGammaSetMarkerTGraphAsym(graphPionRejection[ifilltrue][icalo][ieta], markerStyleEta[ieta], markerSizeEta[ieta], colorEta[ieta] , colorEta[ieta]);
            graphPionRejection[ifilltrue][icalo][ieta]->Draw("samelpe");
          }
          cout << " done l:" << __LINE__ << endl;
        } else if(partEtaCalo[ieta] >= etarangemin[icalo] && partEtaCalo[ieta+1] <= etarangemax[icalo]){
          cout << partEtaCalo[ieta] << " < eta  < " << partEtaCalo[ieta+1];
          if(histPionRejection[ifilltrue][imatch][icalo][ieta]){
            graphPionRejection[ifilltrue][icalo][ieta] = new TGraphAsymmErrors(histPionRejection[ifilltrue][imatch][icalo][ieta]);
            while(graphPionRejection[ifilltrue][icalo][ieta]->GetN()>0 && graphPionRejection[ifilltrue][icalo][ieta]->GetY()[graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1] == -1) graphPionRejection[ifilltrue][icalo][ieta]->RemovePoint(graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1);
            while(graphPionRejection[ifilltrue][icalo][ieta]->GetN()>0 && graphPionRejection[ifilltrue][icalo][ieta]->GetX()[graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1] > energymax[icalo]) graphPionRejection[ifilltrue][icalo][ieta]->RemovePoint(graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1);
            DrawGammaSetMarkerTGraphAsym(graphPionRejection[ifilltrue][icalo][ieta], markerStyleEta[ieta], markerSizeEta[ieta], colorEta[ieta] , colorEta[ieta]);
            graphPionRejection[ifilltrue][icalo][ieta]->Draw("samelpe");
          }
          cout << " done l:" << __LINE__ << endl;
        }
      }
      TLegend* legendPi0Merge           = GetAndSetLegend2(0.65, 0.92, 0.85, 0.92-(5*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
      for (Int_t ieta = 0; ieta < nEta+1; ieta++){
        if(ieta==nEta){
          cout << nominalEtaRegion[icalo][0] << " < eta2  < " << nominalEtaRegion[icalo][1] << endl;
          if(graphPionRejection[ifilltrue][icalo][ieta]){
              legendPi0Merge->AddEntry(graphPionRejection[ifilltrue][icalo][ieta],Form("%1.1f < #eta < %1.1f", nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1]),"p");
          }
        } else if(partEtaCalo[ieta] >= etarangemin[icalo] && partEtaCalo[ieta+1] <= etarangemax[icalo]){
          cout << partEtaCalo[ieta] << " < eta2  < " << partEtaCalo[ieta+1] << endl;
          if(graphPionRejection[ifilltrue][icalo][ieta]){
              legendPi0Merge->AddEntry(graphPionRejection[ifilltrue][icalo][ieta],Form("%1.1f < #eta < %1.1f", partEtaCalo[ieta],partEtaCalo[ieta+1]),"p");
          }
        }
      }
      legendPi0Merge->Draw();
      float toplegval = 0.90;
      drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
      drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
      drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
      // drawLatexAdd("#pi^{#pm}  rej. = (#it{E}/#it{p})_{#pi}<1.6#sigma_{E} / (#it{E}/#it{p})_{#pi} > 0",0.15,toplegval-0.15,textSizeLabelsRel,false,false,false);
      // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

      // DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


      canvasMergingFraction->Update();
      canvasMergingFraction->Print(Form("%s/EoP_pionrejection_%s_%s.%s",outputDir.Data(),caloName[icalo].Data(),str_EoverP_FillTrue[ifilltrue].Data(),suffix.Data()));
    }

    histo2DAccEff->DrawCopy();
    TLegend* legendPi0Merge           = GetAndSetLegend2(0.15, 0.92-0.15, 0.3, 0.92-0.15-(5*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    for (Int_t icalo = 0; icalo < maxcalo; icalo++){
      cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
      if(histPionRejection[ifilltrue][imatch][icalo][nEta]){
        DrawGammaSetMarkerTGraphAsym(graphPionRejection[ifilltrue][icalo][nEta], markerStylePID[icalo], markerSizePID[icalo], colorPID[icalo] , colorPID[icalo]);
        graphPionRejection[ifilltrue][icalo][nEta]->Draw("samecx0");
      }
      cout << " done l:" << __LINE__ << endl;
      legendPi0Merge->AddEntry(graphPionRejection[ifilltrue][icalo][nEta],Form("%s",caloName[icalo].Data()),"p");
    }
    legendPi0Merge->Draw();
    float toplegval = 0.90;
    drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
    drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
    // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

    // DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


    canvasMergingFraction->Update();
    canvasMergingFraction->Print(Form("%s/EoP_pionrejection_AllCalo_%s.%s",outputDir.Data(),str_EoverP_FillTrue[ifilltrue].Data(),suffix.Data()));
  }
    histo2DAccEff->GetYaxis()->SetRangeUser(1,99999);
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{track} (GeV/#it{c})", "#pi^{#pm}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
  histo2DAccEff->DrawCopy();
  Style_t markerStyleCalo[maxcalo] = {20, 47, 33};
  Style_t markerStyleCalo2[maxcalo] = {24, 46, 27};
  Style_t markerStyleCalo3[maxcalo] = {25, 42, 28};
  Size_t markerSizeCalo[maxcalo] = {1.5, 1.8, 2.2};
  TLegend* legendPi0Merge1           = GetAndSetLegend2(0.65, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  TLegend* legendPi0Merge2           = GetAndSetLegend2(0.75, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  TLegend* legendPi0Merge3           = GetAndSetLegend2(0.80, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histPionRejection[1][imatch][icalo][nEta]){
      DrawGammaSetMarkerTGraphAsym(graphPionRejection[1][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
      graphPionRejection[1][icalo][nEta]->SetLineStyle(icalo+1);
      graphPionRejection[1][icalo][nEta]->SetLineWidth(3);
      graphPionRejection[1][icalo][nEta]->Draw("samecx0");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge1->AddEntry(graphPionRejection[1][icalo][nEta]," ","l");
    // legendPi0Merge->AddEntry(graphPionRejection[1][icalo][nEta],Form("%s",caloName[icalo].Data()),"p");
  }
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histPionRejection[0][imatch][icalo][nEta]){
      DrawGammaSetMarkerTGraphAsym(graphPionRejection[0][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
      graphPionRejection[0][icalo][nEta]->Draw("samepe");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge2->AddEntry(graphPionRejection[0][icalo][nEta]," ","p");
    legendPi0Merge3->AddEntry((TObject*)0,Form("%s",caloName[icalo].Data()),"");
  }
  legendPi0Merge1->Draw();
  legendPi0Merge2->Draw();
  legendPi0Merge3->Draw();
  float toplegval = 0.90;
  drawLatexAdd("#it{p}_{rec}",0.75,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#it{p}_{true}",0.65,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  // DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


  canvasMergingFraction->Update();
  canvasMergingFraction->Print(Form("%s/EoP_pionrejection_AllCalo_truerec.%s",outputDir.Data(),suffix.Data()));
  histo2DAccEff->DrawCopy();
  double energymax_elec[maxcalo] = {30, 20, 30.0};

  legendPi0Merge1           = GetAndSetLegend2(0.65, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge2           = GetAndSetLegend2(0.75, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge3           = GetAndSetLegend2(0.80, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histPionRejection[1][imatch][icalo][nEta]){
      DrawGammaSetMarkerTGraphAsym(graphPionRejection[1][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
      graphPionRejection[1][icalo][nEta]->SetLineStyle(icalo+1);
      graphPionRejection[1][icalo][nEta]->SetLineWidth(3);
      graphPionRejection[1][icalo][nEta]->Draw("samecx0");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge1->AddEntry(graphPionRejection[1][icalo][nEta]," ","l");
    // legendPi0Merge->AddEntry(graphPionRejection[1][icalo][nEta],Form("%s",caloName[icalo].Data()),"p");
  }
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histPionRejection_95[1][imatch][icalo][nEta]){
      graphPionRejection_95[1][icalo][nEta] = new TGraphAsymmErrors(histPionRejection_95[1][imatch][icalo][nEta]);
      while(graphPionRejection_95[1][icalo][nEta]->GetN()>0 && graphPionRejection_95[1][icalo][nEta]->GetY()[graphPionRejection_95[1][icalo][nEta]->GetN()-1] == -1) graphPionRejection_95[1][icalo][nEta]->RemovePoint(graphPionRejection_95[1][icalo][nEta]->GetN()-1);
      while(graphPionRejection_95[1][icalo][nEta]->GetN()>0 && graphPionRejection_95[1][icalo][nEta]->GetX()[graphPionRejection_95[1][icalo][nEta]->GetN()-1] > energymax_elec[icalo]) graphPionRejection_95[1][icalo][nEta]->RemovePoint(graphPionRejection_95[1][icalo][nEta]->GetN()-1);
      DrawGammaSetMarkerTGraphAsym(graphPionRejection_95[1][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID_light[icalo] , colorPID_light[icalo]);
      graphPionRejection_95[1][icalo][nEta]->SetLineStyle(icalo+1);
      graphPionRejection_95[1][icalo][nEta]->SetLineWidth(3);
      graphPionRejection_95[1][icalo][nEta]->Draw("samepeX0");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge2->AddEntry(graphPionRejection_95[1][icalo][nEta]," ","p");
    legendPi0Merge3->AddEntry((TObject*)0,Form("%s",caloName[icalo].Data()),"");
  }
  legendPi0Merge1->Draw();
  legendPi0Merge2->Draw();
  legendPi0Merge3->Draw();
  toplegval = 0.90;
  drawLatexAdd("#sigma",0.67,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#epsilon_{95%}",0.75,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  // DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


  canvasMergingFraction->Update();
  canvasMergingFraction->Print(Form("%s/EoP_pionrejection_AllCalo_95comp.%s",outputDir.Data(),suffix.Data()));
  histo2DAccEff->DrawCopy();

  legendPi0Merge1           = GetAndSetLegend2(0.65, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge2           = GetAndSetLegend2(0.75, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge3           = GetAndSetLegend2(0.80, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histPionRejection[1][imatch][icalo][nEta]){
      DrawGammaSetMarkerTGraphAsym(graphPionRejection[1][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
      graphPionRejection[1][icalo][nEta]->SetLineStyle(icalo+1);
      graphPionRejection[1][icalo][nEta]->SetLineWidth(3);
      graphPionRejection[1][icalo][nEta]->Draw("samecx0");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge1->AddEntry(graphPionRejection[1][icalo][nEta]," ","l");
    // legendPi0Merge->AddEntry(graphPionRejection[1][icalo][nEta],Form("%s",caloName[icalo].Data()),"p");
  }
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histPionRejection_fwhm[1][imatch][icalo][nEta]){
      graphPionRejection_fwhm[1][icalo][nEta] = new TGraphAsymmErrors(histPionRejection_fwhm[1][imatch][icalo][nEta]);
      while(graphPionRejection_fwhm[1][icalo][nEta]->GetN()>0 && graphPionRejection_fwhm[1][icalo][nEta]->GetY()[graphPionRejection_fwhm[1][icalo][nEta]->GetN()-1] == -1) graphPionRejection_fwhm[1][icalo][nEta]->RemovePoint(graphPionRejection_fwhm[1][icalo][nEta]->GetN()-1);
      while(graphPionRejection_fwhm[1][icalo][nEta]->GetN()>0 && graphPionRejection_fwhm[1][icalo][nEta]->GetX()[graphPionRejection_fwhm[1][icalo][nEta]->GetN()-1] > energymax_elec[icalo]) graphPionRejection_fwhm[1][icalo][nEta]->RemovePoint(graphPionRejection_fwhm[1][icalo][nEta]->GetN()-1);
      DrawGammaSetMarkerTGraphAsym(graphPionRejection_fwhm[1][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID_light[icalo] , colorPID_light[icalo]);
      graphPionRejection_fwhm[1][icalo][nEta]->SetLineStyle(icalo+1);
      graphPionRejection_fwhm[1][icalo][nEta]->SetLineWidth(3);
      graphPionRejection_fwhm[1][icalo][nEta]->Draw("samepeX0");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge2->AddEntry(graphPionRejection_fwhm[1][icalo][nEta]," ","p");
    legendPi0Merge3->AddEntry((TObject*)0,Form("%s",caloName[icalo].Data()),"");
  }
  legendPi0Merge1->Draw();
  legendPi0Merge2->Draw();
  legendPi0Merge3->Draw();
  toplegval = 0.90;
  drawLatexAdd("#sigma",0.67,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#epsilon_{95%}",0.75,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  // DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


  canvasMergingFraction->Update();
  canvasMergingFraction->Print(Form("%s/EoP_pionrejection_AllCalo_fwhmcomp.%s",outputDir.Data(),suffix.Data()));



  TCanvas* canvaseleceffi       = new TCanvas("canvaseleceffi", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvaseleceffi,  0.1, 0.01, 0.015, 0.095);
  // canvaseleceffi->SetLogy(1);
  canvaseleceffi->SetLogx(1);
  TGraphAsymmErrors* graphElectronEfficiency[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
  TGraphAsymmErrors* graphElectronEfficiency_fwhm[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
  TH2F * histoDummyEffi;
  histoDummyEffi                = new TH2F("histoDummyEffi", "histoDummyEffi",1000, 0.41, 34.9, 1000, 0, 1.29 );
  SetStyleHistoTH2ForGraphs( histoDummyEffi, "#it{p}_{track} (GeV/#it{c})", "electron survival probability",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
    // histoDummyEffi->GetYaxis()->SetRangeUser(1,30999);
  histoDummyEffi->DrawCopy();
  legendPi0Merge1           = GetAndSetLegend2(0.65, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge2           = GetAndSetLegend2(0.75, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge3           = GetAndSetLegend2(0.80, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histElectronEfficiency[1][imatch][icalo][nEta]){
      graphElectronEfficiency[1][icalo][nEta] = new TGraphAsymmErrors(histElectronEfficiency[1][imatch][icalo][nEta]);
      while(graphElectronEfficiency[1][icalo][nEta]->GetN()>0 && graphElectronEfficiency[1][icalo][nEta]->GetY()[graphElectronEfficiency[1][icalo][nEta]->GetN()-1] == -1) graphElectronEfficiency[1][icalo][nEta]->RemovePoint(graphElectronEfficiency[1][icalo][nEta]->GetN()-1);
      while(graphElectronEfficiency[1][icalo][nEta]->GetN()>0 && graphElectronEfficiency[1][icalo][nEta]->GetX()[graphElectronEfficiency[1][icalo][nEta]->GetN()-1] > energymax_elec[icalo]) graphElectronEfficiency[1][icalo][nEta]->RemovePoint(graphElectronEfficiency[1][icalo][nEta]->GetN()-1);
      DrawGammaSetMarkerTGraphAsym(graphElectronEfficiency[1][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
      graphElectronEfficiency[1][icalo][nEta]->SetLineStyle(icalo+1);
      graphElectronEfficiency[1][icalo][nEta]->SetLineWidth(3);
      graphElectronEfficiency[1][icalo][nEta]->Draw("samepeX0");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge1->AddEntry(graphElectronEfficiency[1][icalo][nEta]," ","p");
    // legendPi0Merge->AddEntry(graphElectronEfficiency[1][icalo][nEta],Form("%s",caloName[icalo].Data()),"p");
  }
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histElectronEfficiency_fwhm[1][imatch][icalo][nEta]){
      graphElectronEfficiency_fwhm[1][icalo][nEta] = new TGraphAsymmErrors(histElectronEfficiency_fwhm[1][imatch][icalo][nEta]);
      while(graphElectronEfficiency_fwhm[1][icalo][nEta]->GetN()>0 && graphElectronEfficiency_fwhm[1][icalo][nEta]->GetY()[graphElectronEfficiency_fwhm[1][icalo][nEta]->GetN()-1] == -1) graphElectronEfficiency_fwhm[1][icalo][nEta]->RemovePoint(graphElectronEfficiency_fwhm[1][icalo][nEta]->GetN()-1);
      while(graphElectronEfficiency_fwhm[1][icalo][nEta]->GetN()>0 && graphElectronEfficiency_fwhm[1][icalo][nEta]->GetX()[graphElectronEfficiency_fwhm[1][icalo][nEta]->GetN()-1] > energymax_elec[icalo]) graphElectronEfficiency_fwhm[1][icalo][nEta]->RemovePoint(graphElectronEfficiency_fwhm[1][icalo][nEta]->GetN()-1);
      DrawGammaSetMarkerTGraphAsym(graphElectronEfficiency_fwhm[1][icalo][nEta], markerStyleCalo2[icalo], 1.5*markerSizeCalo[icalo], colorPID_light[icalo] , colorPID_light[icalo]);
      graphElectronEfficiency_fwhm[1][icalo][nEta]->SetLineStyle(icalo+1);
      graphElectronEfficiency_fwhm[1][icalo][nEta]->SetLineWidth(3);
      graphElectronEfficiency_fwhm[1][icalo][nEta]->Draw("samepeX0");
    }
    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge2->AddEntry(graphElectronEfficiency_fwhm[1][icalo][nEta]," ","p");
    legendPi0Merge3->AddEntry((TObject*)0,Form("%s",caloName[icalo].Data()),"");
  }
  legendPi0Merge1->Draw();
  legendPi0Merge2->Draw();
  legendPi0Merge3->Draw();
  drawLatexAdd("FWHM",0.72,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#sigma",0.67,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  DrawGammaLines(0.41, 34.9, 1,1,2,kGray+1, 2);
  DrawGammaLines(0.41, 34.9, 0.95,0.95,2,kGray, 5);


  canvaseleceffi->Update();
  canvaseleceffi->Print(Form("%s/EoP_electroneffi_AllCalo_sigmafwhm.%s",outputDir.Data(),suffix.Data()));


  TCanvas* canvasEPcutValue       = new TCanvas("canvasEPcutValue", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasEPcutValue,  0.1, 0.01, 0.015, 0.095);
  // canvasEPcutValue->SetLogy(1);
  canvasEPcutValue->SetLogx(1);
  TH2F * histoEPcutValue;
  histoEPcutValue                = new TH2F("histoEPcutValue", "histoEPcutValue",1000, 0.51, 49, 1000, 0.5, 1.19 );
  SetStyleHistoTH2ForGraphs( histoEPcutValue, "#it{p}_{track,true} (GeV/#it{c})", "#it{E}/#it{p} cut value",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
  histoEPcutValue->GetYaxis()->SetLabelOffset(0.001);
  histoEPcutValue->GetXaxis()->SetNoExponent();
  histoEPcutValue->GetXaxis()->SetMoreLogLabels(kTRUE);
  TGraphErrors* graphEPcutValue[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
  TGraphErrors* graphEPcutValue_fwhm[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
  TGraphErrors* graphEPcutValue_95[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};

  histoEPcutValue->DrawCopy();
  legendPi0Merge1           = GetAndSetLegend2(0.50, 0.14, 0.70, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge2           = GetAndSetLegend2(0.60, 0.14, 0.80, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge3           = GetAndSetLegend2(0.70, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  TLegend* legendPi0Merge4           = GetAndSetLegend2(0.75, 0.14, 0.95, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    int iEta = nEta;
    int ifilltrue = 1;
    graphEPcutValue[1][icalo][iEta]   = new TGraphErrors(15,partE,EP_cut[ifilltrue][imatch][icalo][nEta]);
    while(graphEPcutValue[1][icalo][iEta]->GetN()>0 && graphEPcutValue[1][icalo][iEta]->GetX()[graphEPcutValue[1][icalo][iEta]->GetN()-1] >= energymax_elec[icalo]) graphEPcutValue[1][icalo][iEta]->RemovePoint(graphEPcutValue[1][icalo][iEta]->GetN()-1);
    DrawGammaSetMarkerTGraph(graphEPcutValue[1][icalo][iEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
    graphEPcutValue[1][icalo][iEta]->SetLineStyle(icalo+1);
    graphEPcutValue[1][icalo][iEta]->SetLineWidth(3);
    graphEPcutValue[1][icalo][iEta]->Draw("samepe");
    
    graphEPcutValue_fwhm[1][icalo][iEta]   = new TGraphErrors(15,partE,EP_cut_FWHM[ifilltrue][imatch][icalo][nEta]);
    while(graphEPcutValue_fwhm[1][icalo][iEta]->GetN()>0 && graphEPcutValue_fwhm[1][icalo][iEta]->GetX()[graphEPcutValue_fwhm[1][icalo][iEta]->GetN()-1] >= energymax_elec[icalo]) graphEPcutValue_fwhm[1][icalo][iEta]->RemovePoint(graphEPcutValue_fwhm[1][icalo][iEta]->GetN()-1);
    DrawGammaSetMarkerTGraph(graphEPcutValue_fwhm[1][icalo][iEta], markerStyleCalo2[icalo], 1.5*markerSizeCalo[icalo], colorPID_light[icalo] , colorPID_light[icalo]);
    graphEPcutValue_fwhm[1][icalo][iEta]->SetLineStyle(icalo+1);
    graphEPcutValue_fwhm[1][icalo][iEta]->SetLineWidth(3);
    graphEPcutValue_fwhm[1][icalo][iEta]->Draw("samepe");
    
    graphEPcutValue_95[1][icalo][iEta]   = new TGraphErrors(15,partE,EP_cut_95[ifilltrue][imatch][icalo][nEta]);
    while(graphEPcutValue_95[1][icalo][iEta]->GetN()>0 && graphEPcutValue_95[1][icalo][iEta]->GetX()[graphEPcutValue_95[1][icalo][iEta]->GetN()-1] >= energymax_elec[icalo]) graphEPcutValue_95[1][icalo][iEta]->RemovePoint(graphEPcutValue_95[1][icalo][iEta]->GetN()-1);
    DrawGammaSetMarkerTGraph(graphEPcutValue_95[1][icalo][iEta], markerStyleCalo3[icalo], 1.5*markerSizeCalo[icalo], colorPID_light2[icalo] , colorPID_light2[icalo]);
    graphEPcutValue_95[1][icalo][iEta]->SetLineStyle(icalo+1);
    graphEPcutValue_95[1][icalo][iEta]->SetLineWidth(3);
    graphEPcutValue_95[1][icalo][iEta]->Draw("samecx0");

    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge1->AddEntry(graphEPcutValue[1][icalo][iEta]," ","p");
    legendPi0Merge2->AddEntry(graphEPcutValue_fwhm[1][icalo][iEta]," ","p");
    legendPi0Merge3->AddEntry(graphEPcutValue_95[1][icalo][iEta]," ","l");
    legendPi0Merge4->AddEntry((TObject*)0,Form("%s",caloName[icalo].Data()),"");
  }

  legendPi0Merge1->Draw();
  legendPi0Merge2->Draw();
  legendPi0Merge3->Draw();
  legendPi0Merge4->Draw();
  // float toplegval = 0.90;
  // drawLatexAdd("#it{p}_{rec}",0.75,0.30,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#it{p}_{true}",0.65,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("1.6#cdot",0.55,0.35,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#sigma",0.52,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#scale[0.75]{FWHM}",0.58,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#epsilon_{95%}",0.70,0.30,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  DrawGammaLines(0.51, 49, 1,1,2,kGray+1, 2);


  canvasEPcutValue->Update();
  canvasEPcutValue->Print(Form("%s/EoP_cutVlue_AllCalo_true.%s",outputDir.Data(),suffix.Data()));


  histoEPcutValue->DrawCopy();
  legendPi0Merge1           = GetAndSetLegend2(0.50, 0.14, 0.70, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge2           = GetAndSetLegend2(0.60, 0.14, 0.80, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge3           = GetAndSetLegend2(0.70, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  legendPi0Merge4           = GetAndSetLegend2(0.75, 0.14, 0.95, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    int iEta = nEta;
    graphEPcutValue[0][icalo][iEta]   = new TGraphErrors(15,partE,EP_cut[0][imatch][icalo][nEta]);
    while(graphEPcutValue[0][icalo][iEta]->GetN()>0 && graphEPcutValue[0][icalo][iEta]->GetX()[graphEPcutValue[0][icalo][iEta]->GetN()-1] >= energymax_elec[icalo]) graphEPcutValue[0][icalo][iEta]->RemovePoint(graphEPcutValue[0][icalo][iEta]->GetN()-1);
    DrawGammaSetMarkerTGraph(graphEPcutValue[0][icalo][iEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
    graphEPcutValue[0][icalo][iEta]->SetLineStyle(icalo+1);
    graphEPcutValue[0][icalo][iEta]->SetLineWidth(3);
    graphEPcutValue[0][icalo][iEta]->Draw("samepe");
    
    graphEPcutValue_fwhm[0][icalo][iEta]   = new TGraphErrors(15,partE,EP_cut_FWHM[0][imatch][icalo][nEta]);
    while(graphEPcutValue_fwhm[0][icalo][iEta]->GetN()>0 && graphEPcutValue_fwhm[0][icalo][iEta]->GetX()[graphEPcutValue_fwhm[0][icalo][iEta]->GetN()-1] >= energymax_elec[icalo]) graphEPcutValue_fwhm[0][icalo][iEta]->RemovePoint(graphEPcutValue_fwhm[0][icalo][iEta]->GetN()-1);
    DrawGammaSetMarkerTGraph(graphEPcutValue_fwhm[0][icalo][iEta], markerStyleCalo2[icalo], 1.5*markerSizeCalo[icalo], colorPID_light[icalo] , colorPID_light[icalo]);
    graphEPcutValue_fwhm[0][icalo][iEta]->SetLineStyle(icalo+1);
    graphEPcutValue_fwhm[0][icalo][iEta]->SetLineWidth(3);
    graphEPcutValue_fwhm[0][icalo][iEta]->Draw("samepe");
    
    graphEPcutValue_95[0][icalo][iEta]   = new TGraphErrors(15,partE,EP_cut_95[0][imatch][icalo][nEta]);
    while(graphEPcutValue_95[0][icalo][iEta]->GetN()>0 && graphEPcutValue_95[0][icalo][iEta]->GetX()[graphEPcutValue_95[0][icalo][iEta]->GetN()-1] >= energymax_elec[icalo]) graphEPcutValue_95[0][icalo][iEta]->RemovePoint(graphEPcutValue_95[0][icalo][iEta]->GetN()-1);
    DrawGammaSetMarkerTGraph(graphEPcutValue_95[0][icalo][iEta], markerStyleCalo3[icalo], 1.5*markerSizeCalo[icalo], colorPID_light2[icalo] , colorPID_light2[icalo]);
    graphEPcutValue_95[0][icalo][iEta]->SetLineStyle(icalo+1);
    graphEPcutValue_95[0][icalo][iEta]->SetLineWidth(3);
    graphEPcutValue_95[0][icalo][iEta]->Draw("samecx0");

    cout << " done l:" << __LINE__ << endl;
    legendPi0Merge1->AddEntry(graphEPcutValue[0][icalo][iEta]," ","p");
    legendPi0Merge2->AddEntry(graphEPcutValue_fwhm[0][icalo][iEta]," ","p");
    legendPi0Merge3->AddEntry(graphEPcutValue_95[0][icalo][iEta]," ","l");
    legendPi0Merge4->AddEntry((TObject*)0,Form("%s",caloName[icalo].Data()),"");
  }

  legendPi0Merge1->Draw();
  legendPi0Merge2->Draw();
  legendPi0Merge3->Draw();
  legendPi0Merge4->Draw();
  // float toplegval = 0.90;
  // drawLatexAdd("#it{p}_{rec}",0.75,0.30,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#it{p}_{true}",0.65,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("1.6#cdot",0.55,0.35,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#sigma",0.52,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#scale[0.75]{FWHM}",0.58,0.30,textSizeLabelsRel,false,false,false);
  drawLatexAdd("#epsilon_{95%}",0.70,0.30,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  DrawGammaLines(0.51, 49, 1,1,2,kGray+1, 2);


  canvasEPcutValue->Update();
  canvasEPcutValue->Print(Form("%s/EoP_cutVlue_AllCalo_rec.%s",outputDir.Data(),suffix.Data()));



  TCanvas* canvasExampleBin       = new TCanvas("canvasExampleBin", "", 200, 10, 1300, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasExampleBin,  0.105, 0.01, 0.015, 0.09);
  // canvasExampleBin->SetLogy(1);
  // canvasExampleBin->SetLogx(1);
  cout << __LINE__ << endl;
  textSizeLabelsRel        = 58./1200;
  TH2F * histoExampleBin;
  histoExampleBin                = new TH2F("histoExampleBin", "histoExampleBin",1000, 0.0, 1.23, 1000, -0.002, 0.19 );
  SetStyleHistoTH2ForGraphs( histoExampleBin, "#it{E}/#it{p}", "probability",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1, 510, 505);//(#times #epsilon_{pur})
  histoExampleBin->GetYaxis()->SetLabelOffset(0.001);
  histoExampleBin->GetXaxis()->SetNoExponent();
  histoExampleBin->GetXaxis()->SetMoreLogLabels(kTRUE);
  int exampleBinCalo[maxcalo] = {3,4,5};
    Color_t colorPi = kRed+2;
    Color_t colorElec = kBlue+2;
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){

    histoExampleBin->DrawCopy();
    TLegend* legendExmplPlot           = GetAndSetLegend2(0.13, 0.885-(2*textSizeLabelsRel), 0.55, 0.885,0.9*textSizeLabelsPixel,2);
    // TLegend* legendPi0Merge2           = GetAndSetLegend2(0.75, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    // legendPi0Merge3           = GetAndSetLegend2(0.70, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    int ieta = nEta;
  cout << __LINE__ << endl;

    // histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][exampleBinCalo[icalo]] 
    // histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][exampleBinCalo[icalo]]
    TH1D* hPionEoPplot = (TH1D*)histEoP_proj_pions[1][imatch][icalo][ieta][exampleBinCalo[icalo]]->Clone(Form("hPionEoPplot%s",caloName[icalo].Data()));
    TH1D* hElectronEoPplot = (TH1D*)histEoP_proj_electrons[1][imatch][icalo][ieta][exampleBinCalo[icalo]]->Clone(Form("hElectronEoPplot%s",caloName[icalo].Data()));
    hPionEoPplot->Scale(1/hPionEoPplot->GetEntries());
    hElectronEoPplot->Scale(1/hElectronEoPplot->GetEntries());
    DrawGammaSetMarker(hPionEoPplot, 33, 4, colorPi , colorPi);
    hPionEoPplot->Draw("samepe");
    DrawGammaSetMarker(hElectronEoPplot, 47, 3.5, colorElec , colorElec);
    hElectronEoPplot->Draw("samepe");
  cout << __LINE__ << endl;


    TH1D* hPionEoPplot_rec = (TH1D*)histEoP_proj_pions[0][imatch][icalo][ieta][exampleBinCalo[icalo]]->Clone(Form("hPionEoPplot_rec%s",caloName[icalo].Data()));
    TH1D* hElectronEoPplot_rec = (TH1D*)histEoP_proj_electrons[0][imatch][icalo][ieta][exampleBinCalo[icalo]]->Clone(Form("hElectronEoPplot_rec%s",caloName[icalo].Data()));
    hPionEoPplot_rec->Scale(1/hPionEoPplot_rec->GetEntries());
    hElectronEoPplot_rec->Scale(1/hElectronEoPplot_rec->GetEntries());
    DrawGammaSetMarker(hPionEoPplot_rec, 27, 4, colorPi-6 , colorPi-6);
    hPionEoPplot_rec->Draw("samepe");
    DrawGammaSetMarker(hElectronEoPplot_rec, 46, 3.5, colorElec-6 , colorElec-6);
    hElectronEoPplot_rec->Draw("samepe");
  cout << __LINE__ << endl;

    cout << " done l:" << __LINE__ << endl;
    legendExmplPlot->AddEntry(hPionEoPplot,"#pi^{-} (true #it{p})","p");
    legendExmplPlot->AddEntry(hElectronEoPplot,"e^{-} (true #it{p})","p");
    legendExmplPlot->AddEntry(hPionEoPplot_rec,"#pi^{-} (rec. #it{p})","p");
    legendExmplPlot->AddEntry(hElectronEoPplot_rec,"e^{-} (rec. #it{p})","p");
    // legendPi0Merge3->AddEntry(graphPionRejection[1][icalo][iEta],Form("%s",caloName[icalo].Data()),"p");

    legendExmplPlot->Draw();
    // float toplegval = 0.90;
    // drawLatexAdd("#it{p}_{rec}",0.75,0.30,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#it{p}_{true}",0.65,0.30,textSizeLabelsRel,false,false,false);
    drawLatexAdd(perfLabel.Data(),0.945,toplegval,textSizeLabelsRel,false,false,true);
    // drawLatexAdd(collisionSystem.Data(),0.94,toplegval-0.05,textSizeLabelsRel,false,false,true);
    // drawLatexAdd(collisionSystem.Data(),0.94,toplegval-0.05,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%s, w/ mat. infront",caloName[icalo].Data()),0.945,toplegval-0.05,textSizeLabelsRel,false,false,true);
    drawLatexAdd(Form("%1.1f< #it{#eta}< %1.1f, #it{E} = %1.1f GeV",nominalEtaRegion[icalo][0],nominalEtaRegion[icalo][1], (partE[exampleBinCalo[icalo]]+partE[exampleBinCalo[icalo]+1])/2),0.14,toplegval,textSizeLabelsRel,false,false,false);

    // DrawGammaLines(EP_cut[0][imatch][icalo][ieta][exampleBinCalo[icalo]], EP_cut[0][imatch][icalo][nEta][exampleBinCalo[icalo]], 0,0.1,2,kGray+1, 3);
    DrawGammaLines(EP_cut[1][imatch][icalo][ieta][exampleBinCalo[icalo]], EP_cut[1][imatch][icalo][nEta][exampleBinCalo[icalo]], 0,0.1,2,kGray+2, 2);


    canvasExampleBin->Update();
    canvasExampleBin->Print(Form("%s/EoP_ExampleBin_%s.%s",outputDir.Data(),caloName[icalo].Data(),suffix.Data()));
  }



  TFile* fileOutput = new TFile(Form("%s/output_EoverP.root",outputDir.Data()),"RECREATE");
  for(int icalo=0;icalo<maxcalo;icalo++){
    // if (!caloEnabled[icalo]) continue;

    fileOutput->mkdir(Form("%s",caloName[icalo].Data()));
    fileOutput->cd(Form("%s",caloName[icalo].Data()));
    if(graphEPcutValue[0][icalo][nEta])graphEPcutValue[0][icalo][nEta]->Write(Form("graphEPcutValue_recP_sigma_%s",caloName[icalo].Data()));
    if(graphEPcutValue_fwhm[0][icalo][nEta])graphEPcutValue_fwhm[0][icalo][nEta]->Write(Form("graphEPcutValue_recP_fwhm_%s",caloName[icalo].Data()));
    if(graphEPcutValue_95[0][icalo][nEta])graphEPcutValue_95[0][icalo][nEta]->Write(Form("graphEPcutValue_recP_95effi_%s",caloName[icalo].Data()));
    if(graphEPcutValue[1][icalo][nEta])graphEPcutValue[1][icalo][nEta]->Write(Form("graphEPcutValue_trueP_sigma_%s",caloName[icalo].Data()));
    if(graphEPcutValue_fwhm[1][icalo][nEta])graphEPcutValue_fwhm[1][icalo][nEta]->Write(Form("graphEPcutValue_trueP_fwhm_%s",caloName[icalo].Data()));
    if(graphEPcutValue_95[1][icalo][nEta])graphEPcutValue_95[1][icalo][nEta]->Write(Form("graphEPcutValue_trueP_95effi_%s",caloName[icalo].Data()));

    if(graphPionRejection[0][icalo][nEta])graphPionRejection[0][icalo][nEta]->Write(Form("graphPionRejection_recP_sigma_%s",caloName[icalo].Data()));
    if(graphPionRejection[1][icalo][nEta])graphPionRejection[1][icalo][nEta]->Write(Form("graphPionRejection_trueP_sigma_%s",caloName[icalo].Data()));
    if(graphPionRejection_fwhm[0][icalo][nEta])graphPionRejection_fwhm[0][icalo][nEta]->Write(Form("graphPionRejection_recP_fwhm_%s",caloName[icalo].Data()));
    if(graphPionRejection_fwhm[1][icalo][nEta])graphPionRejection_fwhm[1][icalo][nEta]->Write(Form("graphPionRejection_trueP_fwhm_%s",caloName[icalo].Data()));
    if(graphPionRejection_95[0][icalo][nEta])graphPionRejection_95[0][icalo][nEta]->Write(Form("graphPionRejection_recP_95effi_%s",caloName[icalo].Data()));
    if(graphPionRejection_95[1][icalo][nEta])graphPionRejection_95[1][icalo][nEta]->Write(Form("graphPionRejection_trueP_95effi_%s",caloName[icalo].Data()));
    if(graphElectronEfficiency[1][icalo][nEta])graphElectronEfficiency[1][icalo][nEta]->Write(Form("graphElectronEfficiency_trueP_sigma_%s",caloName[icalo].Data()));
    if(graphElectronEfficiency_fwhm[1][icalo][nEta])graphElectronEfficiency_fwhm[1][icalo][nEta]->Write(Form("graphElectronEfficiency_trueP_fwhm_%s",caloName[icalo].Data()));
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
