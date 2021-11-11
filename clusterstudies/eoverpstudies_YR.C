#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 

float getEPcutvalue(int icalo, float partE){
  float EP_cut_strength = 1.6;
  float EP_cut = 1;
  
  if(icalo==0) EP_cut = 1. - EP_cut_strength * std::sqrt( 0.5 * 0.5 + 2.2 * 2.2 / partE ) / 100.;
  if(icalo==1) EP_cut = 1. - EP_cut_strength *std::sqrt( 0.8 * 0.8 + 1.7 * 1.7 / partE ) / 100.;
  if(icalo==2) EP_cut = 1. - EP_cut_strength * std::sqrt( 0.1 * 0.1 + 7.3 * 7.3 / partE ) / 100.;

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
  double eopcutvaluemax[maxcalo] = {2.0, 2.0, 2.0};
  double etarangemin[maxcalo] = {-4, -1.7, 1.3};
  double etarangemax[maxcalo] = {-1.7, 1.3, 4.0};
  double energymax[maxcalo] = {50, 20, 50.0};
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem           = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString perfLabel                 = "ECCE G4 simulation";
  
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
  TH1F* histIntCounts_electrons_EoP[Nusefilltrue][maxmatch][maxcalo][nEta+1]= {{{{NULL}}}};
  TH1F* histIntCounts_pions_all[Nusefilltrue][maxmatch][maxcalo][nEta+1]    = {{{{NULL}}}};
  TH1F* histIntCounts_electrons_all[Nusefilltrue][maxmatch][maxcalo][nEta+1]= {{{{NULL}}}};
  TH1F* histPionRejection[Nusefilltrue][maxmatch][maxcalo][nEta+1]          = {{{{NULL}}}};

  double EP_cut[Nusefilltrue][maxmatch][maxcalo][nEta+1][50] = {0};
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
        histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_electrons_EoP_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_pions_all_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histIntCounts_electrons_all_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);
        histPionRejection[ifilltrue][imatch][icalo][ieta] = new TH1F(Form("histPionRejection_%s_%f",caloName[icalo].Data(),minEtaProj),"",nEne-1,partE);

        for(int ebin=0;ebin<nEne-1;ebin++){

          h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->SetRange(  h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin]),h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin+1]));
          h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->SetRange(  h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(minEtaProj),h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(maxEtaProj));
          h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->SetRange( h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin]),h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetYaxis()->FindBin(partE[ebin+1]));
          h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->SetRange( h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(minEtaProj),h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->GetZaxis()->FindBin(maxEtaProj));

          histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("%1.1f_%1.1f_%d_%s_x",minEtaProj,partE[ebin],imatch,caloName[icalo].Data()));
          histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]    = (TH1F*)h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("%1.1f_%1.1f_%d_%s_x",minEtaProj,partE[ebin],imatch,caloName[icalo].Data()));
          // histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_pions_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("histEoP_proj_pions_%s_%f_%f",caloName[icalo].Data(),minEtaProj, partE[ebin]));
          // histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_electrons_EoverP_E_Eta[imatch][icalo][ifilltrue]->Project3D(Form("histEoP_proj_pions_%s_%f_%f",caloName[icalo].Data(),minEtaProj, partE[ebin]));
          float EP_cut_tmp = getEPcutvalue(icalo,(partE[ebin]+partE[ebin+1])/2);
          double EP_shift =0;
          if(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->GetEntries()>10){
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]     = new TF1(Form("fit_%1.1f_%1.1f_%d_%s_x",minEtaProj,partE[ebin],imatch,caloName[icalo].Data()), "gaus", 0.8, 1.05);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParLimits(1, 0.8, 1.0);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParameter(1, 0.95);
            double fitsigmatemp = 0.02;
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParLimits(2, 0.6*fitsigmatemp,3*fitsigmatemp);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParameter(2, fitsigmatemp);
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->Fit(EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec],"L0RMEQ","",0.8, 1.05);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParLimits(1, EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1)*0.9, EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1)*1.1);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParameter(1, EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1));
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParLimits(2, EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(2)*0.5,EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(2)*1.5);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetParameter(2, EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(2));
            histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->Fit(EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec],"L0RMEQ","",EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1)-EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(2)*2, EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1)+EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(2)*3);
            EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->SetRange(EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1)-EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(2)*2, EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1)+EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(2)*3);
            EP_shift     = 1 - EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParameter(1);
            // meanErr  = EP_fit[ifilltrue][imatch][icalo][ieta][padebinelec]->GetParError(1);
            // EP_shift = 1-histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->GetMaximum();
            EP_cut[ifilltrue][imatch][icalo][ieta][padebinelec] = EP_cut_tmp-EP_shift;
          }else {
            // if(padebinelec==0)
              EP_cut[ifilltrue][imatch][icalo][ieta][padebinelec] = EP_cut_tmp;
            // else
            //   EP_cut[ifilltrue][imatch][icalo][ieta][padebinelec] = EP_cut[ifilltrue][imatch][icalo][ieta][padebinelec-1];
          }
          // cout << "EPcutvalue: " << EP_cut[ifilltrue][imatch][icalo][ieta][padebin] << " = " << EP_cut_tmp << " - " << EP_shift << endl;

          Double_t intpions  = histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->Integral(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(EP_cut[ifilltrue][imatch][icalo][ieta][padebin]), histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          if(intpions==0)intpions=1;
          histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpions);

          Double_t intpionsall  = histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->Integral(histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(0.0), histEoP_proj_pions[ifilltrue][imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(1.9));
          histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpionsall);
          padebin++;


          Double_t intelectrons  = histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->Integral(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(EP_cut[ifilltrue][imatch][icalo][ieta][padebin]), histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectrons);

          Double_t intelectronsall  = histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->Integral(histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(0.0), histEoP_proj_electrons[ifilltrue][imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(1.9));
          histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta]->SetBinContent(histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectronsall);

          padebinelec++;

          // cout << (partE[ebin]+partE[ebin+1])/2 << "\tpions: " << histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << "\telectrons: " << histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_electrons_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << endl;
          // cout << "\tallpions: " << histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << "\tallelectrons: " << histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_electrons_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << endl;
          double rejection = histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) / histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2));
          if(isnan(rejection)) rejection = -1;
          if(isinf(rejection)) rejection = -1;
          // cout << (partE[ebin]+partE[ebin+1])/2 << "\trejection: " << rejection << "\tpions: " << histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_EoP[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << "\tallpions: " << histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_all[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2))<< endl;
          histPionRejection[ifilltrue][imatch][icalo][ieta]->SetBinContent(histPionRejection[ifilltrue][imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),rejection);
        }

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
  histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.51, 49, 1000, 1, 59999 );
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{track} (GeV/#it{c})", "#pi^{#pm}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
  histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
  histo2DAccEff->GetXaxis()->SetNoExponent();
  histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
  TGraphAsymmErrors* graphPionRejection[Nusefilltrue][maxcalo][nEta+1]={{{NULL}}};
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
          cout << " done" << endl;
        } else if(partEtaCalo[ieta] >= etarangemin[icalo] && partEtaCalo[ieta+1] <= etarangemax[icalo]){
          cout << partEtaCalo[ieta] << " < eta  < " << partEtaCalo[ieta+1];
          if(histPionRejection[ifilltrue][imatch][icalo][ieta]){
            graphPionRejection[ifilltrue][icalo][ieta] = new TGraphAsymmErrors(histPionRejection[ifilltrue][imatch][icalo][ieta]);
            while(graphPionRejection[ifilltrue][icalo][ieta]->GetN()>0 && graphPionRejection[ifilltrue][icalo][ieta]->GetY()[graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1] == -1) graphPionRejection[ifilltrue][icalo][ieta]->RemovePoint(graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1);
            while(graphPionRejection[ifilltrue][icalo][ieta]->GetN()>0 && graphPionRejection[ifilltrue][icalo][ieta]->GetX()[graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1] > energymax[icalo]) graphPionRejection[ifilltrue][icalo][ieta]->RemovePoint(graphPionRejection[ifilltrue][icalo][ieta]->GetN()-1);
            DrawGammaSetMarkerTGraphAsym(graphPionRejection[ifilltrue][icalo][ieta], markerStyleEta[ieta], markerSizeEta[ieta], colorEta[ieta] , colorEta[ieta]);
            graphPionRejection[ifilltrue][icalo][ieta]->Draw("samelpe");
          }
          cout << " done" << endl;
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
      cout << " done" << endl;
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
    histo2DAccEff->GetYaxis()->SetRangeUser(1,10999);
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{track} (GeV/#it{c})", "#pi^{#pm}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
  histo2DAccEff->DrawCopy();
  Style_t markerStyleCalo[maxcalo] = {20, 47, 33};
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
    cout << " done" << endl;
    legendPi0Merge1->AddEntry(graphPionRejection[1][icalo][nEta]," ","l");
    // legendPi0Merge->AddEntry(graphPionRejection[1][icalo][nEta],Form("%s",caloName[icalo].Data()),"p");
  }
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << nominalEtaRegion[icalo][0] << " < eta  < " << nominalEtaRegion[icalo][1];
    if(histPionRejection[0][imatch][icalo][nEta]){
      DrawGammaSetMarkerTGraphAsym(graphPionRejection[0][icalo][nEta], markerStyleCalo[icalo], 1.5*markerSizeCalo[icalo], colorPID[icalo] , colorPID[icalo]);
      graphPionRejection[0][icalo][nEta]->Draw("samepe");
    }
    cout << " done" << endl;
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
}
