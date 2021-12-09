#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>
#include <TDatabasePDG.h>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 


double mysine(double* x, double* par) { double 

Amplitude = par[0];
double wavelength = par[1];
double phase = par[2];
return Amplitude*sin(2*TMath::Pi()/wavelength*x[0]+phase); 

}

TF1* FitExpPlusGaussian(TH1F* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t icalo, Double_t ptcenter, Int_t iDataMC);
TF1* FitExpPlusGaussianPol2(TH1F* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t icalo, Double_t ptcenter, Int_t iDataMC);

double GaussExpLinear(double* x, double* par) {
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double p5 = par[5];
  double p6 = par[6];

  double func;
  if(x[0]<p1){
    func = (p0*(TMath::Exp(-0.5*pow(((x[0]-p1)/p2),2))+TMath::Exp((x[0]-p1)/p3)*(1.-TMath::Exp(-0.5*pow(((x[0]-p1)/p2),2))))+p4+p5*x[0]);

  }else{
    func = (p0*(TMath::Exp(-0.5*pow(((x[0]-p1)/p2),2))+TMath::Exp(-(x[0]-p1)/p6)*(1.-TMath::Exp(-0.5*pow(((x[0]-p1)/p1),2))))+p4+p5*x[0]);
  }
  return func;
}

void pi0resolution(
    TString inputFileName     = "central_allculsters.root",
    TString clusterizerName   =  "MA",
    TString suffix            = "pdf",
    Bool_t doFitting          = kTRUE,
    TString collisionsSys     = "SinglePart"
){  

  const int maxcalo = 3;
  // TString caloName[maxcalo+1] = {"FEMC", "EEMC", "BECAL", "Hybrid"};
  TString caloName[maxcalo+1] = { "EEMC", "BECAL", "FEMC", "Hybrid"};
  TString caloNamePlot[maxcalo+1] = { "EEMC", "BEMC", "FEMC", "Hybrid"};
  bool caloactive[maxcalo] = {false};
  Double_t mass_pi0PDG = 0.1349770;

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem           = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString perfLabel                 = "#it{#bf{ECCE}} simulation";
  
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


  
  TString outputDir                 = Form("plotsResolution");
  if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFile                      = new TFile(inputFileName.Data());
  if(inputFile->IsZombie()){
    cout << "inputFile not found!" << endl;
  }

  //************************** Read data **************************************************
  const int nEne = 26;
  // const int nEta = 40;
  const int nPhi = 25;
  
  const static Double_t partE[]   = {   0.7, 0.9, 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 4.5,  5,
                                        5.5, 6, 7, 8, 10,
                                       12.5,15,17.5,20,22.5,25,  30,35,  40, 50};
  double energymax = 40.;

  TH2F* histEtrueMinv[maxcalo+1] = {NULL};
  TH2F* histEtrueMinv_allcls[maxcalo+1] = {NULL};
  TH2F* histEtrueMinv_allcalo= NULL;
  TH2F* histEtrueMinv_allcalo_unmatched= NULL;
  TH2F* histpTMinv= NULL;
  TH2F* histErecMinv[maxcalo+1] = {NULL};
  TH2F* histErecMinv_allcls[maxcalo+1] = {NULL};
  TH2F* histErecMinv_allcalo = NULL;
  TH2F* histErecMinv_allcalo_unmatched= NULL;

  TH1F* histEtrueMinvbin[maxcalo+1][50]            = {NULL};
  TH1F* histEtrueMinvbin_allcls[maxcalo+1][50]            = {NULL};
  TH1F* histEtrueMinvbin_allcalo[50]            = {NULL};

  TF1*  fithistEtrueMinv[maxcalo+1][50] = {NULL};
  TF1*  fithistEtrueMinv_allclus[50] = {NULL};
  TF1*  fithistEtrueMinv_allcalo[50] = {NULL};
  TF1*  fithistEtrueMinv2[50] = {NULL};
  TF1*  fithistEtrueMinvFINAL[maxcalo+1][50] = {NULL};
  TF1*  fithistEtrueMinvFINAL_allcalo[50] = {NULL};

  // TH1F* histResoEVsMinv = new TH1F("histResoEVsMinv","",100,0,41) = {NULL};
  TH1F* histResoEVsMinv[maxcalo+1] = {NULL};
  // TH1F* histMeanEVsMinv = new TH1F("histMeanEVsMinv","",100,0,41) = {NULL};
  TH1F* histMeanEVsMinv[maxcalo+1] = {NULL};
  // TH1F* histSigmaEVsMinv = new TH1F("histSigmaEVsMinv","",100,0,41) = {NULL};
  TH1F* histSigmaEVsMinv[maxcalo+1] = {NULL};
  TH1F* histResoEVsMinv_allcalo= NULL;
  TH1F* histMeanEVsMinv_allcalo= NULL;
  TH1F* histSigmaEVsMinv_allcalo= NULL;
  
  TH1F* histEtrueMinvbinSlice[50]            = {NULL};
  TF1*  fithistEtrueMinvSlice[50] = {NULL};
  TF1*  fithistEtrueMinvSliceFINAL[50] = {NULL};
  
  TH3F* histEtaEMinv= NULL;
  TH1F* histEEtaMinvbin[50][50]         = {{NULL}};
  TF1* fitEEtaMinvbin[50][50]           = {{NULL}};
  TF1* fitEEtaMinvbinFINAL[50][50]           = {{NULL}};
  TH1F* histResoEEtaVsMinv[50];
  TH1F* histMeanEEtaVsMinv[50];
  TH1F* histSigmaEEtaVsMinv[50];

  TH2F* histEClusters;
  TH1F* histEClusterbin[50] = {NULL};
  TH1F* histEProb[50] = {NULL};
  TH2F* histM02p1Pt= NULL; 
  TH2F* histM02p2Pt= NULL; 
  TH2F* histM02SCPt= NULL;
  TH2F* histM20p1Pt= NULL; 
  TH2F* histM20p2Pt= NULL;
  TH2F* histM20SCPt= NULL;


  TH1F* h_Pi0Spec_separatecls[maxcalo+1] = {NULL};
  TH1F* h_Pi0Spec_mergedclus[maxcalo+1] = {NULL};
  TH1F* h_Pi0Spec_total[maxcalo+1] = {NULL};


  Int_t nParticles                  = 7;
  TString readNames[7]              = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_", ""};
  TString labelPart[7]              = {"#gamma", "e^{#pm}", "#pi^{#pm}", "p/#bar{p}", "K^{#pm}", "n", "all"};
  Color_t colorPart[7]              = {807, kBlue+1, kRed+2, kGreen+2, kViolet+2, kGray+2, kBlack};
  Style_t markerStylePart[7]        = {20, 46, 33, 30, 28, 42, 25};
  Style_t lineStylePart[7]          = {3, 4, 5, 7, 9, 3, 1};
  Size_t markerSizePart[7]          = {2.5*1.5, 2*1.5, 3*1.5, 2.7*1.5, 2*1.5, 3*1.5, 2*1.5};

  TH1F* hist1DResoEtabins3D[100][100] = {{NULL}};
  TH3F* hist_E_eta_Clusters;

  Int_t nbins = 100; 
  double ebin = 1.0;
  int padebin = 0; 
  double minmassfit = 0.05;
  double maxmassfit = 0.15;
  double mean, meanErr, sigma, sigmaErr;
  double sigmat = 2.3548;
  int clusmax = 10;

  double Ptbin = 1;
  double Ptmax = 30;
  double binpT = 0.5;
  double ptskip = 1.;
  TH1F* histPtMinvbin[50];
  TH3F* hist_Pt_eta_Clusters; 
  TH1F* histPtProb[50];
  TH1F* hist1DResoPttabins3D[50][50] = {{NULL}};

  
  // vector<double> allmean;
  // vector<double> allsigma;
  std::array<std::vector<double>, maxcalo> allmean;
  std::array<std::vector<double>, maxcalo> allsigma;
  std::vector<double> allmean_allcalo;
  std::vector<double> allsigma_allcalo;

  int padptbin = 0;
  TF1* fithistPtMinv[50];
  TF1* fithistPtMinvFINAL[50];
  TH1F* histResoPtVsMinv = new TH1F("histResoPtVsMinv","",100,0,31);
  TH1F* histMeanPtVsMinv = new TH1F("histMeanPtVsMinv","",100,0,31);
  TH1F* histSigmaPtVsMinv = new TH1F("histSigmaPtVsMinv","",100,0,31);


  TString NameErecon_Minv_allcalo = Form("h_Erec_Minv_allcalo");
  histErecMinv_allcalo    = (TH2F*)inputFile->Get(NameErecon_Minv_allcalo.Data());
  if(!histErecMinv_allcalo){ cout << "\t" <<NameErecon_Minv_allcalo.Data() << " not found" << endl;}
  TString NameE_Minv_allcalo = Form("h_Etrue_Minv_allcalo" );
  histEtrueMinv_allcalo    = (TH2F*)inputFile->Get(NameE_Minv_allcalo.Data());
  if(!histEtrueMinv_allcalo){ cout << "\t" <<NameE_Minv_allcalo.Data() << " not found" << endl; }

  TString NameErecon_Minv_allcalo_unmatched = Form("h_Erec_Minv_allcalo_unmatched");
  histErecMinv_allcalo_unmatched    = (TH2F*)inputFile->Get(NameErecon_Minv_allcalo_unmatched.Data());
  if(!histErecMinv_allcalo_unmatched){ cout << "\t" <<NameErecon_Minv_allcalo_unmatched.Data() << " not found" << endl;}
  TString NameE_Minv_allcalo_unmatched = Form("h_Etrue_Minv_allcalo_unmatched" );
  histEtrueMinv_allcalo_unmatched    = (TH2F*)inputFile->Get(NameE_Minv_allcalo_unmatched.Data());
  if(!histEtrueMinv_allcalo_unmatched){ cout << "\t" <<NameE_Minv_allcalo_unmatched.Data() << " not found" << endl; }

  bool allcaloplotted = false;
  bool allcalofitted = false;

  for(int icalo=0;icalo<maxcalo;icalo++){
    cout << "processing " << caloName[icalo].Data() << endl;
    TString NameE_Minv = Form("%s/h_Etrue_Minv_allcls_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    // TString NameE_Minv = Form("%s/h_Etrue_Minv_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    histEtrueMinv[icalo]    = (TH2F*)inputFile->Get(NameE_Minv.Data());
    if(!histEtrueMinv[icalo]){
      cout << "\t" <<NameE_Minv.Data() << " not found" << endl;
      continue;
    } else {
      caloactive[icalo]=true;
    };
    TString NameErecon_Minv = Form("%s/h_Erec_Minv_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    histErecMinv[icalo]    = (TH2F*)inputFile->Get(NameErecon_Minv.Data());
    if(!histErecMinv[icalo]){ cout << "\t" <<NameErecon_Minv.Data() << " not found" << endl;}
    // TString NameE_Minv_allcls = Form("%s/h_Etrue_Minv_allcls_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    // histEtrueMinv_allcls[icalo]    = (TH2F*)inputFile->Get(NameE_Minv_allcls.Data());
    // if(!histEtrueMinv_allcls[icalo]){ cout << "\t" <<NameE_Minv_allcls.Data() << " not found" << endl; }
    TString NameErecon_Minv_allcls = Form("%s/h_Erec_Minv_allcls_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    histErecMinv_allcls[icalo]    = (TH2F*)inputFile->Get(NameErecon_Minv_allcls.Data());
    if(!histErecMinv_allcls[icalo]){ cout << "\t" <<NameErecon_Minv_allcls.Data() << " not found" << endl; }

    cout << "\tinputs loaded" << endl;

    if (histEtrueMinv[icalo]->GetEntries() > 100 || histErecMinv[icalo]->GetEntries()>100){
      TCanvas* Minvpi02D = new TCanvas(Form("Minvpi02D%s",caloName[icalo].Data()),"",0,0,1500,800);
      DrawGammaCanvasSettings( Minvpi02D, 0.095, 0.01, 0.01, 0.105);
      Minvpi02D->Divide(2,1);
      Minvpi02D->cd(1);
      histEtrueMinv[icalo]->GetXaxis()->SetRangeUser(0.,40);
      histEtrueMinv[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      // histEtrueMinv_allcls[icalo]->GetXaxis()->SetRangeUser(0.,40);
      // histEtrueMinv_allcls[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      SetStyleHistoTH2ForGraphs(histEtrueMinv[icalo],"E^{MC}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      histEtrueMinv[icalo]->Draw("colz");
      Minvpi02D->cd(2);
      histErecMinv[icalo]->GetXaxis()->SetRangeUser(0.,40);
      histErecMinv[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      SetStyleHistoTH2ForGraphs(histErecMinv[icalo],"E^{reco}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      histErecMinv[icalo]->Draw("colz");
      
      TLatex Tl1;
      Tl1.SetTextSize(0.035);
      Tl1.DrawLatex(0.55*energymax,0.04,perfLabel);
      Tl1.DrawLatex(0.55*energymax,0.03,Form("Calorimeter: %s ", caloNamePlot[icalo].Data()));
      Tl1.DrawLatex(0.55*energymax,0.02,Form("Clusterization type: %s ", clusterizerName.Data()));
      Tl1.DrawLatex(0.55*energymax,0.01,"#pi ^{0} #rightarrow #gamma #gamma");
      Tl1.Draw();
      
      Minvpi02D->Print(Form("%s/Mass_vs_E_%s_%s.%s", outputDir.Data(), caloNamePlot[icalo].Data(), clusterizerName.Data(), suffix.Data()));

      cout << "\tMass_vs_E plot saved" << endl;
      if(!allcaloplotted){
        Minvpi02D = new TCanvas(Form("Minvpi02Dallcalo"),"",0,0,1500,800);
        DrawGammaCanvasSettings( Minvpi02D, 0.095, 0.01, 0.01, 0.105);
        Minvpi02D->Divide(2,1);
        Minvpi02D->cd(1);
        histEtrueMinv_allcalo->GetXaxis()->SetRangeUser(0.,40);
        histEtrueMinv_allcalo->GetYaxis()->SetRangeUser(0.,0.2);
        SetStyleHistoTH2ForGraphs(histEtrueMinv_allcalo,"E^{MC}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
        histEtrueMinv_allcalo->Draw("colz");
        Minvpi02D->cd(2);
        histErecMinv_allcalo->GetXaxis()->SetRangeUser(0.,40);
        histErecMinv_allcalo->GetYaxis()->SetRangeUser(0.,0.2);
        SetStyleHistoTH2ForGraphs(histErecMinv_allcalo,"E^{reco}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
        histErecMinv_allcalo->Draw("colz");
        
        TLatex Tl1;
        Tl1.SetTextSize(0.035);
        Tl1.DrawLatex(0.55*energymax,0.04,perfLabel);
        Tl1.DrawLatex(0.55*energymax,0.03,Form("Hybrid Method (all calos)"));
        Tl1.DrawLatex(0.55*energymax,0.02,Form("Clusterization type: %s ", clusterizerName.Data()));
        Tl1.DrawLatex(0.55*energymax,0.01,"#pi ^{0} #rightarrow #gamma #gamma");
        Tl1.Draw();
        
        Minvpi02D->Print(Form("%s/Mass_vs_E_%s_%s.%s", outputDir.Data(), "Hybrid", clusterizerName.Data(), suffix.Data()));
        allcaloplotted = true;
        cout << "\tMass_vs_E plot for allcalo saved" << endl;
  
        Minvpi02D = new TCanvas(Form("Minvpi02Dallcalo_unmatched"),"",0,0,1500,800);
        DrawGammaCanvasSettings( Minvpi02D, 0.095, 0.01, 0.01, 0.105);
        Minvpi02D->Divide(2,1);
        Minvpi02D->cd(1);
        histEtrueMinv_allcalo_unmatched->GetXaxis()->SetRangeUser(0.,40);
        histEtrueMinv_allcalo_unmatched->GetYaxis()->SetRangeUser(0.,0.2);
        SetStyleHistoTH2ForGraphs(histEtrueMinv_allcalo_unmatched,"E^{MC}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
        histEtrueMinv_allcalo_unmatched->Draw("colz");
        Minvpi02D->cd(2);
        histErecMinv_allcalo_unmatched->GetXaxis()->SetRangeUser(0.,40);
        histErecMinv_allcalo_unmatched->GetYaxis()->SetRangeUser(0.,0.2);
        SetStyleHistoTH2ForGraphs(histErecMinv_allcalo_unmatched,"E^{reco}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
        histErecMinv_allcalo_unmatched->Draw("colz");
        
        Tl1.SetTextSize(0.035);
        Tl1.DrawLatex(0.55*energymax,0.04,perfLabel);
        Tl1.DrawLatex(0.55*energymax,0.03,Form("Hybrid Method (all calos)"));
        Tl1.DrawLatex(0.55*energymax,0.02,Form("Clusterization type: %s ", clusterizerName.Data()));
        Tl1.DrawLatex(0.55*energymax,0.01,"#pi ^{0} #rightarrow #gamma #gamma");
        Tl1.Draw();
        
        Minvpi02D->Print(Form("%s/Mass_vs_E_%s_%s.%s", outputDir.Data(), "Hybrid_unmatched", clusterizerName.Data(), suffix.Data()));
        allcaloplotted = true;
        cout << "\tMass_vs_E plot for allcalo unmatched saved" << endl;
      }
      // while(ebin < energymax){
      padebin = 0;

      histResoEVsMinv[icalo] = new TH1F(Form("histResoEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
      histMeanEVsMinv[icalo] = new TH1F(Form("histMeanEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
      histSigmaEVsMinv[icalo] = new TH1F(Form("histSigmaEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
      for(int ebin=0;ebin<nEne-1;ebin++){
        histEtrueMinvbin[icalo][padebin]    = (TH1F*)histEtrueMinv[icalo]->ProjectionY(Form("projection_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]), 
                                              histEtrueMinv[icalo]->GetXaxis()->FindBin(partE[ebin]), histEtrueMinv[icalo]->GetXaxis()->FindBin(partE[ebin+1]),"e");
        // histEtrueMinvbin_allcls[icalo][padebin]    = (TH1F*)histEtrueMinv_allcls[icalo]->ProjectionY(Form("projection_Minv_allcls_%s_%f",caloName[icalo].Data(), partE[ebin]), 
        //                                       histEtrueMinv_allcls[icalo]->GetXaxis()->FindBin(partE[ebin]), histEtrueMinv_allcls[icalo]->GetXaxis()->FindBin(partE[ebin+1]),"e");

        if(histEtrueMinvbin[icalo][padebin]->GetEntries() > 10){
          if(icalo==2){
            histEtrueMinvbin[icalo][padebin]->Rebin(6);
            // histEtrueMinvbin_allcls[icalo][padebin]->Rebin(6);
          } else {
            histEtrueMinvbin[icalo][padebin]->Rebin(4);
            // histEtrueMinvbin_allcls[icalo][padebin]->Rebin(4);
          }
          mean = meanErr = sigma = sigmaErr = 0;

          if (doFitting){   

            fithistEtrueMinv[icalo][padebin] = FitExpPlusGaussian (histEtrueMinvbin[icalo][padebin], minmassfit, maxmassfit, icalo, (partE[ebin]+partE[ebin+1])/2, 0);
            // fithistEtrueMinv[icalo][padebin] = FitExpPlusGaussianPol2 (histEtrueMinvbin[icalo][padebin], minmassfit, maxmassfit, icalo, (partE[ebin]+partE[ebin+1])/2, 0);
            // fithistEtrueMinv[icalo][padebin]    = new TF1(Form("fit_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]), "crystalball", 0.05, 0.15);
            // fithistEtrueMinv[icalo][padebin]->SetParameters(0.55*histEtrueMinvbin[icalo][padebin]->GetMaximum(),histEtrueMinvbin[icalo][padebin]->GetMean(),2.*histEtrueMinvbin[icalo][padebin]->GetRMS(),2.*histEtrueMinvbin[icalo][padebin]->GetRMS(),2.5*histEtrueMinvbin[icalo][padebin]->GetRMS());
            histEtrueMinvbin[icalo][padebin]->Fit(fithistEtrueMinv[icalo][padebin],"L0RMEQ","",minmassfit, maxmassfit);

            // fithistEtrueMinvFINAL[icalo][padebin] = histEtrueMinvbin[icalo][padebin]->GetFunction(Form("fit_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]));
            fithistEtrueMinvFINAL[icalo][padebin] = fithistEtrueMinv[icalo][padebin];

            mean     = fithistEtrueMinvFINAL[icalo][padebin]->GetParameter(1);
            meanErr  = fithistEtrueMinvFINAL[icalo][padebin]->GetParError(1);
            sigma    = fithistEtrueMinvFINAL[icalo][padebin]->GetParameter(2);
            sigmaErr = fithistEtrueMinvFINAL[icalo][padebin]->GetParError(2);
            allmean[icalo].push_back(mean);
            allsigma[icalo].push_back(sigma);
          }


          Double_t reso     = sigma/mean;
          Double_t resoErr  = TMath::Sqrt( TMath::Power(sigmaErr/sigma,2) + TMath::Power(meanErr/mean,2) );

          histResoEVsMinv[icalo]->SetBinContent(histResoEVsMinv[icalo]->FindBin(partE[ebin]),  100*reso);
          histResoEVsMinv[icalo]->SetBinError(histResoEVsMinv[icalo]->FindBin(partE[ebin]),  100*resoErr);  
          histMeanEVsMinv[icalo]->SetBinContent(histMeanEVsMinv[icalo]->FindBin(partE[ebin]), mean);
          histMeanEVsMinv[icalo]->SetBinError(histMeanEVsMinv[icalo]->FindBin(partE[ebin]), meanErr);

          histSigmaEVsMinv[icalo]->SetBinContent(histSigmaEVsMinv[icalo]->FindBin(partE[ebin]), sigma/sigmat);
          histSigmaEVsMinv[icalo]->SetBinError(histSigmaEVsMinv[icalo]->FindBin(partE[ebin]), sigmaErr/sigmat);
          padebin++;
        }
        // ebin = ebin + 2.; 
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
          Tl.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", caloNamePlot[icalo].Data()));
          Tl.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
          Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEtrueMinvbin[icalo][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
          DrawGammaSetMarker(histEtrueMinvbin[icalo][padnumber],34,1, kBlack, kBlack);
          histEtrueMinvbin[icalo][padnumber]->GetYaxis()->SetNdivisions(5);
          histEtrueMinvbin[icalo][padnumber]->GetXaxis()->SetNdivisions(5);
          histEtrueMinvbin[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
          histEtrueMinvbin[icalo][padnumber]->Draw("ep");
          // DrawGammaSetMarker(histEtrueMinvbin_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
          // histEtrueMinvbin_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
          // histEtrueMinvbin_allcls[icalo][padnumber]->Draw("same,hist");
          if(fithistEtrueMinvFINAL[icalo][padnumber]){
            fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
            fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
            fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
            DrawGammaLines(allmean[icalo][padnumber], allmean[icalo][padnumber], 0, 0.5*fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(),2, 2);
            DrawGammaLines(allmean[icalo][padnumber] + allsigma[icalo][padnumber], allmean[icalo][padnumber]+ allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
            DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
            if(z==0){
              legend1->AddEntry(histEtrueMinvbin[icalo][padnumber], "Simulation", "p");
              legend1->AddEntry(fithistEtrueMinvFINAL[icalo][padnumber], "Fit", "l");
            }
        }
          drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
          padnumber++;
        }
      }

      mass2D->Print(Form("%s/Mass_%s_%s.%s", outputDir.Data(), caloNamePlot[icalo].Data(), clusterizerName.Data(), suffix.Data()));


      if(!allcalofitted){
        // while(ebin < energymax){
        int padebin = 0;

        histResoEVsMinv_allcalo = new TH1F(Form("histResoEVsMinv_%s","hybrid"),"",nEne-1,partE);
        histMeanEVsMinv_allcalo = new TH1F(Form("histMeanEVsMinv_%s","hybrid"),"",nEne-1,partE);
        histSigmaEVsMinv_allcalo = new TH1F(Form("histSigmaEVsMinv_%s","hybrid"),"",nEne-1,partE);
        for(int ebin=0;ebin<nEne-1;ebin++){
          histEtrueMinvbin_allcalo[padebin]    = (TH1F*)histEtrueMinv_allcalo->ProjectionY(Form("projection_Minv_%s_%f","hybrid", partE[ebin]), 
                                                histEtrueMinv_allcalo->GetXaxis()->FindBin(partE[ebin]), histEtrueMinv_allcalo->GetXaxis()->FindBin(partE[ebin+1]),"e");

          if(histEtrueMinvbin_allcalo[padebin]->GetEntries() > 10){
            histEtrueMinvbin_allcalo[padebin]->Rebin(4);
            mean = meanErr = sigma = sigmaErr = 0;

            if (doFitting){

              // fithistEtrueMinv_allcalo[padebin]    = new TF1(Form("fit_Minv_%s_%f","hybrid", partE[ebin]), "crystalball", 0.05, 0.15);
              // fithistEtrueMinv_allcalo[padebin]->SetParameters(0.55*histEtrueMinvbin_allcalo[padebin]->GetMaximum(),histEtrueMinvbin_allcalo[padebin]->GetMean(),2.*histEtrueMinvbin_allcalo[padebin]->GetRMS(),2.*histEtrueMinvbin_allcalo[padebin]->GetRMS(),2.5*histEtrueMinvbin_allcalo[padebin]->GetRMS());
              fithistEtrueMinv_allcalo[padebin] = FitExpPlusGaussian (histEtrueMinvbin_allcalo[padebin], minmassfit, maxmassfit, icalo, (partE[ebin]+partE[ebin+1])/2, 0);
              // fithistEtrueMinv_allcalo[padebin] = FitExpPlusGaussianPol2 (histEtrueMinvbin_allcalo[padebin], minmassfit, maxmassfit, icalo, (partE[ebin]+partE[ebin+1])/2, 0);
              histEtrueMinvbin_allcalo[padebin]->Fit(fithistEtrueMinv_allcalo[padebin],"L0RMEQ","",minmassfit, maxmassfit);
              
              fithistEtrueMinvFINAL_allcalo[padebin] = fithistEtrueMinv_allcalo[padebin];
              // fithistEtrueMinvFINAL_allcalo[padebin] = histEtrueMinvbin_allcalo[padebin]->GetFunction(Form("fit_Minv_%s_%f","hybrid", partE[ebin]));
                  
              mean     = fithistEtrueMinvFINAL_allcalo[padebin]->GetParameter(1);
              meanErr  = fithistEtrueMinvFINAL_allcalo[padebin]->GetParError(1);
              sigma    = fithistEtrueMinvFINAL_allcalo[padebin]->GetParameter(2);
              sigmaErr = fithistEtrueMinvFINAL_allcalo[padebin]->GetParError(2);
              allmean_allcalo.push_back(mean);
              allsigma_allcalo.push_back(sigma);
            }


            Double_t reso     = sigma/mean;
            Double_t resoErr  = TMath::Sqrt( TMath::Power(sigmaErr/sigma,2) + TMath::Power(meanErr/mean,2) );
              
            histResoEVsMinv_allcalo->SetBinContent(histResoEVsMinv_allcalo->FindBin(partE[ebin]),  100*reso);
            histResoEVsMinv_allcalo->SetBinError(histResoEVsMinv_allcalo->FindBin(partE[ebin]),  100*resoErr);  
            histMeanEVsMinv_allcalo->SetBinContent(histMeanEVsMinv_allcalo->FindBin(partE[ebin]), mean);
            histMeanEVsMinv_allcalo->SetBinError(histMeanEVsMinv_allcalo->FindBin(partE[ebin]), meanErr);

            histSigmaEVsMinv_allcalo->SetBinContent(histSigmaEVsMinv_allcalo->FindBin(partE[ebin]), sigma/sigmat);
            histSigmaEVsMinv_allcalo->SetBinError(histSigmaEVsMinv_allcalo->FindBin(partE[ebin]), sigmaErr/sigmat);
            padebin++;
          }
          // ebin = ebin + 2.; 
        }


        mass2D = new TCanvas("mass2D","",0,0,1000,800);
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetTitleBorderSize(0);
        gStyle->SetTitleFontSize(0.07); 
        DrawGammaCanvasSettings( mass2D, 0.095, 0.01, 0.01, 0.10);
        mass2D->Divide(5,6);
        int padnumber = 0;
        legend1 = new TLegend(0., 0.01, 0.4, 0.25);
        for(int z=0; z< padebin +1 ; z++ ){
          if(z==4){
            mass2D->cd(z+1);
            TLatex Tl;
            Tl.SetTextSize(0.1);
            Tl.DrawLatex(0.03,0.9,perfLabel);
            Tl.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", "Hybrid"));
            Tl.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
            Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

            legend1->Draw();
          } else{
            mass2D->cd(z+1);
            //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
            SetStyleHistoTH1ForGraphs( histEtrueMinvbin_allcalo[padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
            DrawGammaSetMarker(histEtrueMinvbin_allcalo[padnumber],34,1, kBlack, kBlack);
            histEtrueMinvbin_allcalo[padnumber]->GetYaxis()->SetNdivisions(5);
            histEtrueMinvbin_allcalo[padnumber]->GetXaxis()->SetNdivisions(5);
            histEtrueMinvbin_allcalo[padnumber]->GetXaxis()->SetRangeUser(0,0.3);
            histEtrueMinvbin_allcalo[padnumber]->Draw("ep");
            // DrawGammaSetMarker(histEtrueMinvbin_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
            // histEtrueMinvbin_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
            // histEtrueMinvbin_allcls[icalo][padnumber]->Draw("same,hist");
            fithistEtrueMinvFINAL_allcalo[padnumber]->SetLineColor(kBlue);
            fithistEtrueMinvFINAL_allcalo[padnumber]->SetLineWidth(1);
            fithistEtrueMinvFINAL_allcalo[padnumber]->Draw("same");
            drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
            DrawGammaLines(allmean_allcalo[padnumber], allmean_allcalo[padnumber], 0, 0.5*fithistEtrueMinvFINAL_allcalo[padnumber]->GetMaximum(),2, 2);
            DrawGammaLines(allmean_allcalo[padnumber] + allsigma_allcalo[padnumber], allmean_allcalo[padnumber]+ allsigma_allcalo[padnumber], 0, fithistEtrueMinvFINAL_allcalo[padnumber]->GetMaximum(), 1, 12);
            DrawGammaLines(allmean_allcalo[padnumber] - allsigma_allcalo[padnumber], allmean_allcalo[padnumber]- allsigma_allcalo[padnumber], 0, fithistEtrueMinvFINAL_allcalo[padnumber]->GetMaximum(), 1, 12);
            if(z==0){
              legend1->AddEntry(histEtrueMinvbin_allcalo[padnumber], "Simulation", "p");
              legend1->AddEntry(fithistEtrueMinvFINAL_allcalo[padnumber], "Fit", "l");
            }
            padnumber++;
          }
        }
        mass2D->Print(Form("%s/Mass_%s_%s.%s", outputDir.Data(), "hybrid", clusterizerName.Data(), suffix.Data()));

        allcalofitted = true;
      }

      //============= Fit Results======================================
      TCanvas* fitresults = new TCanvas("fitresults","",0,0,1000,500);
      DrawGammaCanvasSettings( fitresults, 0.15, 0.15, 0.15, 0.15);
    
      fitresults->Divide(2,1);
      fitresults->cd(1);

      histMeanEVsMinv[icalo]->GetYaxis()->SetRangeUser(0.08,0.16);
      SetStyleHistoTH1ForGraphs( histMeanEVsMinv[icalo], "(M_{#gamma #gamma}", "", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.,1.5);
      DrawGammaSetMarker(histMeanEVsMinv[icalo],34,1, kBlack, kBlack);   
      histMeanEVsMinv[icalo]->GetYaxis()->SetRangeUser(0,1.3*histMeanEVsMinv[icalo]->GetMaximum()); 
      histMeanEVsMinv[icalo]->Draw("e,p");
      DrawGammaLines(0, energymax, mass_pi0PDG, mass_pi0PDG,2, 13);

      TLatex Tl3;
      Tl3.SetTextSize(0.035);
      Tl3.DrawLatex(2,0.18,perfLabel);
      Tl3.DrawLatex(2,0.17,Form("Calorimeter: %s ", caloNamePlot[icalo].Data()));
      Tl3.DrawLatex(2,0.16,Form("Clusterization type: %s ", clusterizerName.Data()));
      Tl3.DrawLatex(2,0.14,"#pi ^{0} #rightarrow #gamma #gamma");
      Tl3.Draw();
      fitresults->cd(2);

      histSigmaEVsMinv[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      SetStyleHistoTH1ForGraphs(histSigmaEVsMinv[icalo], "E^{MC}", "Width (M_{#gamma #gamma})", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.3);  
      DrawGammaSetMarker(histSigmaEVsMinv[icalo],34,1, kBlack, kBlack);   
      histSigmaEVsMinv[icalo]->GetYaxis()->SetRangeUser(0,0.015); 
      histSigmaEVsMinv[icalo]->Draw("e,p");
      
      
      fitresults->Print(Form("%s/Mean_Width_%s_%s.%s", outputDir.Data(), caloNamePlot[icalo].Data(), clusterizerName.Data(), suffix.Data()));
    }


    TString Name_mclus = Form("%s/h_Pi0Spec_mergedclus_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    h_Pi0Spec_mergedclus[icalo]    = (TH1F*)inputFile->Get(Name_mclus.Data());
    TString Name_sepclus = Form("%s/h_Pi0Spec_separatecls_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    h_Pi0Spec_separatecls[icalo]    = (TH1F*)inputFile->Get(Name_sepclus.Data());
    TString Name_totclus = Form("%s/h_Pi0Spec_total_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    h_Pi0Spec_total[icalo]    = (TH1F*)inputFile->Get(Name_totclus.Data());
  }

  
  //============= End Fit Results======================================

  float minPtPi0 = 0.6;
  float maxPtPi0woMerged = 49;

  Double_t mesonMassExpectPi0                     = TDatabasePDG::Instance()->GetParticle(111)->Mass();
  Double_t mesonMassExpectEta                     = TDatabasePDG::Instance()->GetParticle(221)->Mass();
  Double_t arrayBoundariesX1_4[2];
  Double_t arrayBoundariesY1_4[3];
  Double_t relativeMarginsX[3];
  Double_t relativeMarginsY[3];
  textSizeLabelsPixel             = 50;
  ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

  TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
  DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

  TPad* padWidthPi0               = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
  DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
  padWidthPi0->Draw();

  TPad* padMassPi0                = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
  DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
  padMassPi0->Draw();

  TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
  DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
  padMassLegend1->SetFillStyle(0);
  padMassLegend1->Draw();

  padWidthPi0->cd();
  padWidthPi0->SetLogx();

  Double_t margin                 = relativeMarginsX[0]*2.7*1350;
  Double_t textsizeLabelsWidth    = 0;
  Double_t textsizeFacWidth       = 0;
  if (padWidthPi0->XtoPixel(padWidthPi0->GetX2()) < padWidthPi0->YtoPixel(padWidthPi0->GetY1())){
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
      textsizeFacWidth            = (Double_t)1./padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
  } else {
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->YtoPixel(padWidthPi0->GetY1());
      textsizeFacWidth            = (Double_t)1./padWidthPi0->YtoPixel(padWidthPi0->GetY1());
  }

  TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, minPtPi0,maxPtPi0woMerged ,1000., -30, 40);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{E} (GeV)", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin),512,505,62,42);//#it{p}_{T} (GeV/#it{c})
  histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(0.1,19.5);//24.5);
  histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
  histo2DAllPi0FWHM->DrawCopy();

  Color_t colorDet[maxcalo+1];
  colorDet[0] = getCaloColor("EEMC", false);
  colorDet[1] = getCaloColor("BEMC", false);
  colorDet[2] = getCaloColor("FEMC", false);
  colorDet[3] = kBlack;
  Style_t markerStyleDet[maxcalo+1]   = {20, 47, 33, 25};//, 28, 30, 42, 46, 24, 25, 27, 28, 30};
  Size_t markerSizeDet[maxcalo+1]     = {2.5, 3, 3 , 2.5};//, 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9,1.5, 1.8 };
  TGraphAsymmErrors* graphPi0FWHMMeV[maxcalo+1];
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    if(caloactive[icalo]){
      if(histSigmaEVsMinv[icalo]){
        histSigmaEVsMinv[icalo]->Scale(1000);
        graphPi0FWHMMeV[icalo] = new TGraphAsymmErrors(histSigmaEVsMinv[icalo]);
        if(graphPi0FWHMMeV[icalo]){
          if(icalo==1){
            while(graphPi0FWHMMeV[icalo]->GetX()[graphPi0FWHMMeV[icalo]->GetN()-1] > 30) graphPi0FWHMMeV[icalo]->RemovePoint(graphPi0FWHMMeV[icalo]->GetN()-1);
          } else if(icalo==2){
            while(graphPi0FWHMMeV[icalo]->GetX()[graphPi0FWHMMeV[icalo]->GetN()-1] > 12) graphPi0FWHMMeV[icalo]->RemovePoint(graphPi0FWHMMeV[icalo]->GetN()-1);
          }
            DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[icalo], markerStyleDet[icalo], markerSizeDet[icalo], colorDet[icalo] , colorDet[icalo]);
            // DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[0] , colorDet[0]);
            graphPi0FWHMMeV[icalo]->Draw("p,same,e");
            // DrawGammaSetMarkerTGraphAsym(graphPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            // graphPi0TrueFWHMMeV[i]->Draw("p,same,e");
        }
      }
    }
  }
  TGraphAsymmErrors* graphPi0FWHMMeV_allcalo;
  if(histSigmaEVsMinv_allcalo){
    histSigmaEVsMinv_allcalo->Scale(1000);
    graphPi0FWHMMeV_allcalo = new TGraphAsymmErrors(histSigmaEVsMinv_allcalo);
    if(graphPi0FWHMMeV_allcalo){
      while(graphPi0FWHMMeV_allcalo->GetX()[graphPi0FWHMMeV_allcalo->GetN()-1] > 30) graphPi0FWHMMeV_allcalo->RemovePoint(graphPi0FWHMMeV_allcalo->GetN()-1);
        DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV_allcalo, markerStyleDet[maxcalo], markerSizeDet[maxcalo], colorDet[maxcalo] , colorDet[maxcalo]);
        // DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[0] , colorDet[0]);
        graphPi0FWHMMeV_allcalo->Draw("p,same,e");
        // DrawGammaSetMarkerTGraphAsym(graphPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
        // graphPi0TrueFWHMMeV[i]->Draw("p,same,e");
    }
  }

  // TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
  // SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
  // labelLegendAMass->SetTextFont(43);
  // labelLegendAMass->Draw();

  TLatex *labelMassPerf       = new TLatex(0.13,0.87,perfLabel.Data());//"ALICE performance");
  SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
  labelMassPerf->SetTextFont(43);
  labelMassPerf->Draw();
  TLatex *labelMassEnergy     = new TLatex(0.13,0.78,Form("%s clusters", clusterizerName.Data()));
  SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4);
  labelMassEnergy->SetTextFont(43);
  labelMassEnergy->Draw();
  TLatex *labelMassPi0        = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4);
  labelMassPi0->SetTextFont(43);
  labelMassPi0->Draw();

  padMassPi0->cd();
  padMassPi0->SetLogx();

  Double_t textsizeLabelsMass         = 0;
  Double_t textsizeFacMass            = 0;
  if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
      textsizeFacMass                 = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
  } else {
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
      textsizeFacMass                 = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
  }

  TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minPtPi0,maxPtPi0woMerged, 1000., 80.1, 167.9);//, 100.1, 160.9);//125.1, 155.9);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{E}_{#pi^{0}} (GeV)", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin),512,505,42,42);
  histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0Mass->GetYaxis()->SetRangeUser(80.1, 149.9);//(131.1, 147.9);//125.1, 155.9);
  histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0Mass->GetXaxis()->SetNoExponent();
  histo2DAllPi0Mass->DrawCopy();
  DrawGammaLines(minPtPi0,maxPtPi0woMerged , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,2, kGray+1,7);
  TGraphAsymmErrors* graphPi0Mass[maxcalo];
  TLegend* legendMassWidth   = GetAndSetLegend2(0.12, 0.2, 0.32, 0.2+4*textsizeLabelsMass, textSizeLabelsPixel);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    if(caloactive[icalo] && histMeanEVsMinv[icalo]){
      histMeanEVsMinv[icalo]->Scale(1000);
      graphPi0Mass[icalo] = new TGraphAsymmErrors(histMeanEVsMinv[icalo]);
      if(icalo==1){
          while(graphPi0Mass[icalo]->GetX()[graphPi0Mass[icalo]->GetN()-1] > 30) graphPi0Mass[icalo]->RemovePoint(graphPi0Mass[icalo]->GetN()-1);
        } else if(icalo==2){
          while(graphPi0Mass[icalo]->GetX()[graphPi0Mass[icalo]->GetN()-1] > 12) graphPi0Mass[icalo]->RemovePoint(graphPi0Mass[icalo]->GetN()-1);
        }
      if(graphPi0Mass[icalo]){
        DrawGammaSetMarkerTGraphAsym(graphPi0Mass[icalo], markerStyleDet[icalo], markerSizeDet[icalo], colorDet[icalo] , colorDet[icalo]);
        graphPi0Mass[icalo]->Draw("p,same,e");
        // DrawGammaSetMarkerTGraphAsym(graphPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
        // graphPi0TrueMass[i]->Draw("p,same,e");
        legendMassWidth->AddEntry(graphPi0Mass[icalo],caloNamePlot[icalo].Data(),"p");
      }
    }
  }
  TGraphAsymmErrors* graphPi0Mass_allcalo;
  if(histMeanEVsMinv_allcalo){
    histMeanEVsMinv_allcalo->Scale(1000);
    graphPi0Mass_allcalo = new TGraphAsymmErrors(histMeanEVsMinv_allcalo);
    if(graphPi0Mass_allcalo){
      while(graphPi0Mass_allcalo->GetX()[graphPi0Mass_allcalo->GetN()-1] > 30) graphPi0Mass_allcalo->RemovePoint(graphPi0Mass_allcalo->GetN()-1);
      DrawGammaSetMarkerTGraphAsym(graphPi0Mass_allcalo, markerStyleDet[maxcalo], markerSizeDet[maxcalo], colorDet[maxcalo] , colorDet[maxcalo]);
      graphPi0Mass_allcalo->Draw("p,same,e");
      // DrawGammaSetMarkerTGraphAsym(graphPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
      // graphPi0TrueMass[i]->Draw("p,same,e");
      legendMassWidth->AddEntry(graphPi0Mass_allcalo,"Hybrid","p");
    }
  }
  legendMassWidth->Draw();


  // TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
  // SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
  // labelLegendBMass->SetTextFont(43);
  // labelLegendBMass->Draw();

  // legendM02Data->SetMargin(0.05/(0.9-0.8));
  // //********************************** Defintion of the Legend **************************************************
  // Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
  // Double_t  rowsLegendMass2[14]= {0.84,0.66,0.50,0.33,0.16,0.00,0.16,0.16,0.16,0.16,0.16,0.16,0.16};
  // //******************* Offsets ***********************
  // Double_t offsetMarkerXMass2         = 0.1;
  // Double_t offsetMarkerYMass2         = 0.1;
  // //****************** Scale factors ******************
  // Double_t scaleMarkerMass2           = 1.2;

  // padMassLegend1->cd();
  // //****************** first Column **************************************************
  // TLatex *textMassPCM[10];
  // Int_t legendRunningIndex = 0;
  // for (Int_t icalo = 0; icalo < maxcalo; icalo++){
  //     if(graphPi0Mass[icalo]){
  //         textMassPCM[icalo]                  = new TLatex(columnsLegendMass2[0],rowsLegendMass2[legendRunningIndex+1],caloName[icalo].Data());
  //         SetStyleTLatex( textMassPCM[icalo], textSizeLabelsPixel,4);
  //         textMassPCM[icalo]->SetTextFont(43);
  //         textMassPCM[icalo]->Draw();
  //         legendRunningIndex++;
  //     }
  // }
  // //****************** second Column *************************************************
  // TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
  // SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
  // textMassData->SetTextFont(43);
  // textMassData->Draw();
  // TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
  // SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
  // textMassMC->SetTextFont(43);
  // textMassMC->Draw();

  // TMarker* markerPCMPi0Mass[10];
  // TMarker* markerPCMPi0MassMC[10];
  // legendRunningIndex = 0;
  // for (Int_t i = 0; i < 13; i++){
  //     if(graphPi0Mass[i] && graphPi0TrueMass[i]){
  //         markerPCMPi0Mass[i]             = CreateMarkerFromGraph(graphPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
  //         markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2);
  //         markerPCMPi0MassMC[i]           = CreateMarkerFromGraph(graphPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
  //         markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2);
  //         legendRunningIndex++;
  //     }
  // }

  canvasMassWidthPi0->Update();
  canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));




    textSizeLabelsPixel                 = 100*3.5/5;
    TCanvas* canvasInvMassSamplePlot    = new TCanvas("canvasInvMassSamplePlot","",0,0,1500,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.105, 0.01, 0.039, 0.093);
    Double_t textSizeLabelPPMass             = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;

    Style_t markerStyleInvMassSGBG      = 0;
    Size_t markerSizeInvMassSGBG        = 0;
    Color_t markerColorInvMassSGBG      = kBlack;
    Style_t markerStyleInvMassMBG       = 24;
    Size_t markerSizeInvMassMBG         = 3;
    Color_t markerColorInvMassMBG       = kGray+2;
    Color_t markerColorInvMassMBG1      = kGray+3;
    Color_t markerColorInvMassMBG2      = kGray+1;
    Style_t markerStyleInvMassBG        = 20;
    Size_t markerSizeInvMassBG          = 3;
    Color_t markerColorInvMassBG        = kBlack;
    Style_t markerStyleInvMassSG        = 20;
    Size_t markerSizeInvMassSG          = 4;
    Color_t markerColorInvMassSG        = kRed+2;
    Color_t fitColorInvMassSG           = kAzure+2;

    Double_t marginInvMass          = 0.1*1500;
    Double_t textsizeLabelsInvMass  = 0;
    Double_t textsizeFacInvMass     = 0;
    if (canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) < canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1())){
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
    } else {
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
    }
    cout << textsizeLabelsInvMass << endl;
    cout << __LINE__ << endl;

    TH2F * histo2DPi0InvMassDummy;
    histo2DPi0InvMassDummy             = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.05,0.21,21000,-1000,200000);
    // histo2DPi0InvMassDummy             = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.05,0.249,21000,-1000,200000);
    SetStyleHistoTH2ForGraphs(histo2DPi0InvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                            0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass),505,510,42,42);
  int exampleBins[maxcalo+1] = {7,7,6,11};
  bool plotinvmass = true;
  if(plotinvmass){
    for (Int_t icalo = 0; icalo < maxcalo+1; icalo++){
      if(histEtrueMinvbin[icalo][exampleBins[icalo]]){
        Int_t iptbinplot = exampleBins[icalo];
        TH1D* histoSignalPlusBGPi0                        = nullptr;
        if(icalo<maxcalo){
          histoSignalPlusBGPi0                        = (TH1D*)histEtrueMinvbin[icalo][iptbinplot]->Clone(Form("Mapping_GG_InvMass_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
        } else {
          histoSignalPlusBGPi0                        = (TH1D*)histEtrueMinvbin_allcalo[iptbinplot]->Clone(Form("Mapping_GG_InvMass_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
        }
        TH1D* histoSignalPi0                              = nullptr;
        TH1D* histoCombBGPi0                              = nullptr;
        TF1* fFitGausExp                                  = nullptr;
        if(icalo<maxcalo){
          fFitGausExp  = (TF1*)fithistEtrueMinvFINAL[icalo][iptbinplot]->Clone(Form("Signal_InvMassFit_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
        } else {
          fFitGausExp = (TF1*)fithistEtrueMinvFINAL_allcalo[iptbinplot]->Clone(Form("Signal_InvMassFit_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
        }
        fFitGausExp->SetNpx(10000);

        canvasInvMassSamplePlot->cd();
        // histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
        histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(1.2*histoSignalPlusBGPi0->GetMinimum(),1.2*histoSignalPlusBGPi0->GetMaximum());
        histo2DPi0InvMassDummy->DrawCopy();

        TLatex *labelInvMassPtRangel = new TLatex(0.945,0.9,Form("#pi^{0}: %2.1f < #it{E}_{#pi^{0}} < %2.1f GeV/#it{c}",partE[iptbinplot],partE[iptbinplot+1]));
        // TLatex *labelInvMassPtRangel = new TLatex(0.945,0.9,Form("#pi^{0}: %2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",partE[iptbinplot],partE[iptbinplot+1]));
        if(histoSignalPlusBGPi0){
          // DrawGammaSetMarker(histoSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
          DrawGammaSetMarker(histoSignalPlusBGPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSGBG, markerColorInvMassSGBG);
          histoSignalPlusBGPi0->SetLineWidth(1);
          // histoSignalPlusBGPi0->Draw("hist,e,same");
          histoSignalPlusBGPi0->Draw("same");
        }
        // if(histoCombBGPi0){
        //   DrawGammaSetMarker(histoCombBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
        //   histoCombBGPi0->SetLineWidth(markerSizeInvMassMBG);
        //   histoCombBGPi0->Draw("same");
        // }
        // if(histoSignalPi0){
        //   DrawGammaSetMarker(histoSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        //   histoSignalPi0->SetStats(kFALSE);
        //   histoSignalPi0->Draw("same");
        // }
        fFitGausExp->SetLineWidth(1);
        fFitGausExp->SetRange(0,0.255);
        fFitGausExp->SetLineColor(fitColorInvMassSG);
        fFitGausExp->Draw("same");

        TLatex *labelInvMassALICE      = new TLatex(0.135,0.9,perfLabel.Data());
        SetStyleTLatex( labelInvMassALICE, 0.85*textSizeLabelsPixel,4);
        labelInvMassALICE->SetTextFont(43);
        labelInvMassALICE->Draw();

        // TLatex *labelInvMassInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textSizeLabelPPMass,collisionSystempp8TeV.Data());
        TLatex *labelInvMassInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.9*textSizeLabelPPMass,collisionSystem.Data());
        SetStyleTLatex( labelInvMassInvMassEnergy, 0.85*textSizeLabelsPixel,4);
        labelInvMassInvMassEnergy->SetTextFont(43);
        labelInvMassInvMassEnergy->Draw();

        TLatex *labelInvMassInvMassTriggerl      = new TLatex(0.135,0.9-0.9*2*0.9*textSizeLabelPPMass,caloName[icalo].Data());
        SetStyleTLatex( labelInvMassInvMassTriggerl, 0.85*textSizeLabelsPixel,4);
        labelInvMassInvMassTriggerl->SetTextFont(43);
        labelInvMassInvMassTriggerl->Draw();

        // TLatex *labelInvMassInvMassReco  = new TLatex(0.135,0.9-0.9*3*0.8*textSizeLabelPPMass,"EMC");
        // SetStyleTLatex( labelInvMassInvMassReco, 0.85*textSizeLabelsPixel,4);
        // labelInvMassInvMassReco->SetTextFont(43);
        // labelInvMassInvMassReco->Draw();

        SetStyleTLatex( labelInvMassPtRangel, 0.85*textSizeLabelsPixel,4);
        labelInvMassPtRangel->SetTextAlign(31);
        labelInvMassPtRangel->SetTextFont(43);
        labelInvMassPtRangel->Draw();

        DrawGammaLines(fFitGausExp->GetParameter(1)-(0.06),fFitGausExp->GetParameter(1)-(0.06),1.2*histoSignalPlusBGPi0->GetMinimum(),0.2*histoSignalPlusBGPi0->GetMaximum(),2,kGray+2,7);
        DrawGammaLines(fFitGausExp->GetParameter(1)+(0.080),fFitGausExp->GetParameter(1)+(0.080),1.2*histoSignalPlusBGPi0->GetMinimum(),0.2*histoSignalPlusBGPi0->GetMaximum(),2,kGray+2,7);
        // TLegend* legendInvMassl  =GetAndSetLegend2(0.62, 0.87-5*0.8*0.75*textSizeLabelPPMass, 0.85, 0.87, 0.85*textSizeLabelsPixel);
        int legendreduction = 0;
        if(!histoCombBGPi0)legendreduction+=2;
        if(!histoSignalPi0)legendreduction++;
        TLegend* legendInvMassl  = GetAndSetLegend2(0.62, 0.87-(5-legendreduction)*0.9*0.75*textSizeLabelPPMass, 0.85, 0.87, 0.85*textSizeLabelsPixel);
        legendInvMassl->SetMargin(0.25);
        legendInvMassl->AddEntry(histoSignalPlusBGPi0,"Raw real events","p");
        // legendInvMassl->AddEntry(histoSignalPlusBGPi0,"Raw real events","l");
        if(histoCombBGPi0){
          legendInvMassl->AddEntry(histoCombBGPi0,"Mixed event +","p");
          legendInvMassl->AddEntry((TObject*)0,"remain. BG","");
        } else {
          // legendInvMassl->AddEntry((TObject*)0,"Mixed event +","");
          // legendInvMassl->AddEntry((TObject*)0,"remain. BG","");
        }
        if(histoSignalPi0){
          legendInvMassl->AddEntry(histoSignalPi0,"BG subtracted","p");
        } else {
          // legendInvMassl->AddEntry((TObject*)0,"BG subtracted","");
        }
        legendInvMassl->AddEntry(fFitGausExp, "Fit","l");
        legendInvMassl->Draw();
        canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBin%d_%s.%s",outputDir.Data(),iptbinplot,caloName[icalo].Data(), suffix.Data()));
      }
    }
  }


    // h_Pi0Spec_mergedclus[icalo]    = (TH2F*)inputFile->Get(Name_mclus.Data());
    // h_Pi0Spec_separatecls[icalo]    = (TH2F*)inputFile->Get(Name_sepclus.Data());
    // h_Pi0Spec_total[icalo]    = (TH2F*)inputFile->Get(Name_totclus.Data());

  TH1F* h_Pi0_mergedfraction[maxcalo] = {NULL};
  TH1F* h_Pi0_separatefraction[maxcalo] = {NULL};

  // **********************************************************************************************************************
  // ******************************** Acceptance * Efficiency for pi0 single measurement 8TeV **************************
  // **********************************************************************************************************************
  textSizeLabelsPixel             = 55;
  textSizeLabelsRel      = 55./1200;
  cout << textSizeLabelsRel << endl;
  cout << __LINE__ << endl;

  TCanvas* canvasMergingFraction       = new TCanvas("canvasMergingFraction", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasMergingFraction,  0.1, 0.01, 0.015, 0.095);
  // canvasMergingFraction->SetLogy(1);
  // canvasMergingFraction->SetLogx(1);
  // float minPtPi0MF = 5;
  TH2F * histo2DAccEff;
  histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, minPtPi0, maxPtPi0woMerged, 1000, 0, 1.09 );
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{E}_{#pi^{0}} (GeV)", "#it{f}_{merg.}",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510,42,42);//(#times #epsilon_{pur})
  histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
  histo2DAccEff->GetXaxis()->SetNoExponent();
  histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAccEff->DrawCopy();
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
      if(h_Pi0Spec_mergedclus[icalo]){
        h_Pi0Spec_mergedclus[icalo]->Sumw2();
        h_Pi0Spec_mergedclus[icalo]->Rebin(4);
        h_Pi0Spec_total[icalo]->Sumw2();
        h_Pi0Spec_total[icalo]->Rebin(4);
        h_Pi0_mergedfraction[icalo] = (TH1F*)h_Pi0Spec_mergedclus[icalo]->Clone(Form("mergedfrachisto_%d",icalo));
        // h_Pi0_mergedfraction[icalo]->Divide(h_Pi0Spec_total[icalo]);
        h_Pi0_mergedfraction[icalo]->Divide(h_Pi0_mergedfraction[icalo],h_Pi0Spec_total[icalo], 1, 1,"B" );
          DrawGammaSetMarker(h_Pi0_mergedfraction[icalo], markerStyleDet[icalo], 1.2*markerSizeDet[icalo], colorDet[icalo] , colorDet[icalo]);
          if(icalo==1)h_Pi0_mergedfraction[icalo]->GetXaxis()->SetRangeUser(2,50);
          if(icalo==0 || icalo==2)h_Pi0_mergedfraction[icalo]->GetXaxis()->SetRangeUser(5,50);
          h_Pi0_mergedfraction[icalo]->Draw("same,p,e");
      }
  }

  TLegend* legendPi0Merge           = GetAndSetLegend2(0.45, 0.13, 0.95, 0.13+(1*textSizeLabelsRel),0.9*textSizeLabelsPixel,3);
  for (Int_t icalo = maxcalo-1; icalo > -1; icalo--){
      if(h_Pi0Spec_mergedclus[icalo]){
          legendPi0Merge->AddEntry(h_Pi0_mergedfraction[icalo],caloNamePlot[icalo].Data(),"p");
      }
  }
  legendPi0Merge->Draw();
  float lowlegval = 0.20;
  drawLatexAdd(perfLabel.Data(),0.95,lowlegval+0.05,textSizeLabelsRel,false,false,true);
  // drawLatexAdd(collisionSystem.Data(),0.95,lowlegval+0.05,textSizeLabelsRel,false,false,true);
  drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.95,lowlegval,textSizeLabelsRel,false,false,true);

  DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


  canvasMergingFraction->Update();
  canvasMergingFraction->Print(Form("%s/Pi0_MergingFraction.%s",outputDir.Data(),suffix.Data()));
 

  cout << __LINE__ << endl;
  TH2F* h_Pi0Spec_mergedclus_Eta[maxcalo] = {NULL};
  TH2F* h_Pi0Spec_separatecls_Eta[maxcalo] = {NULL};
  TH2F* h_Pi0Spec_total_Eta[maxcalo] = {NULL};

  for (Int_t icalo = 0; icalo < maxcalo; icalo++){

    TString Name_mclus_Eta = Form("%s/h_Pi0Spec_mergedclus_Eta_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    h_Pi0Spec_mergedclus_Eta[icalo]    = (TH2F*)inputFile->Get(Name_mclus_Eta.Data());
    TString Name_sepclus_Eta = Form("%s/h_Pi0Spec_separatecls_Eta_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    h_Pi0Spec_separatecls_Eta[icalo]    = (TH2F*)inputFile->Get(Name_sepclus_Eta.Data());
    TString Name_totclus_Eta = Form("%s/h_Pi0Spec_total_Eta_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    h_Pi0Spec_total_Eta[icalo]    = (TH2F*)inputFile->Get(Name_totclus_Eta.Data());
  }
  cout << __LINE__ << endl;

  TH1F* h_Pi0_mergedfraction_Eta[maxcalo][nEta] = {NULL};
  TH1F* h_Pi0_total_Eta[maxcalo][nEta] = {NULL};

  double etarangemin[maxcalo] = {-4, -1.7, 1.3};
  double etarangemax[maxcalo] = {-1.7, 1.3, 4.0};
  float minEtaProj = 0;
  float maxEtaProj = 0;
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << caloName[icalo] << endl;
    for(int ieta=0;ieta<nEta;ieta++){
      // if(ieta>nEta) break;
      minEtaProj = partEtaCalo[ieta];
      maxEtaProj = partEtaCalo[ieta+1];
      cout << ieta << "\t" <<  minEtaProj << " < eta < " << maxEtaProj << endl;

      h_Pi0_mergedfraction_Eta[icalo][ieta]    = (TH1F*)h_Pi0Spec_mergedclus_Eta[icalo]->ProjectionX(Form("h_Pi0_mergedfraction_Eta_%s_%f",caloName[icalo].Data(), minEtaProj), 
                                            h_Pi0Spec_mergedclus_Eta[icalo]->GetYaxis()->FindBin(minEtaProj), h_Pi0Spec_mergedclus_Eta[icalo]->GetYaxis()->FindBin(maxEtaProj),"e");
      h_Pi0_mergedfraction_Eta[icalo][ieta]->Rebin(4);
      h_Pi0_mergedfraction_Eta[icalo][ieta]->Sumw2();
      h_Pi0_total_Eta[icalo][ieta]    = (TH1F*)h_Pi0Spec_total_Eta[icalo]->ProjectionX(Form("h_Pi0_total_Eta_%s_%f",caloName[icalo].Data(), minEtaProj), 
                                            h_Pi0Spec_total_Eta[icalo]->GetYaxis()->FindBin(minEtaProj), h_Pi0Spec_total_Eta[icalo]->GetYaxis()->FindBin(maxEtaProj),"e");
      h_Pi0_total_Eta[icalo][ieta]->Sumw2();
      h_Pi0_total_Eta[icalo][ieta]->Rebin(4);
      h_Pi0_mergedfraction_Eta[icalo][ieta]->Divide(h_Pi0_mergedfraction_Eta[icalo][ieta],h_Pi0_total_Eta[icalo][ieta], 1, 1,"B" );

    }

    histo2DAccEff->DrawCopy();
    TLegend* legendPi0Merge           = GetAndSetLegend2(0.64, 0.13, 0.92, 0.13+(5*0.95*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    for(int ieta=0;ieta<nEta;ieta++){
      if(partEtaCalo[ieta] >= etarangemin[icalo] && partEtaCalo[ieta+1] <= etarangemax[icalo]){
        if(h_Pi0_mergedfraction_Eta[icalo][ieta]){
            DrawGammaSetMarker(h_Pi0_mergedfraction_Eta[icalo][ieta], markerStyleEta[ieta], 1.2*markerSizeEta[ieta], colorEta[ieta] , colorEta[ieta]);
            if(icalo==1)h_Pi0_mergedfraction_Eta[icalo][ieta]->GetXaxis()->SetRangeUser(2,50);
            if(icalo==0 || icalo==2)h_Pi0_mergedfraction_Eta[icalo][ieta]->GetXaxis()->SetRangeUser(5,50);
            h_Pi0_mergedfraction_Eta[icalo][ieta]->Draw("same,p,e");
            legendPi0Merge->AddEntry(h_Pi0_mergedfraction_Eta[icalo][ieta],Form("%1.1f < #it{#eta} < %1.1f",partEtaCalo[ieta],partEtaCalo[ieta+1]),"p");
        }
      }
    }
    legendPi0Merge->Draw();
    float lowlegval = 0.15+(5*0.95*textSizeLabelsRel);
    drawLatexAdd(perfLabel.Data(),0.95,lowlegval+0.05,textSizeLabelsRel,false,false,true);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.95,lowlegval,textSizeLabelsRel,false,false,true);
    // drawLatexAdd(collisionSystem.Data(),0.95,lowlegval+0.05,textSizeLabelsRel,false,false,true);

    DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


    canvasMergingFraction->Update();
    canvasMergingFraction->Print(Form("%s/Pi0_MergingFraction_vsEta_%s.%s",outputDir.Data(),caloName[icalo].Data(),suffix.Data()));
 
  }
}


TF1* FitExpPlusGaussian(TH1F* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t icalo, Double_t ptcenter , Int_t iDataMC){

    Double_t mesonAmplitude = 0;
    for(Int_t i = histo->FindBin(0.07); i < histo->FindBin(0.2) ; i++ ){
      if(histo->GetBinContent(i) > mesonAmplitude) mesonAmplitude = histo->GetBinContent(i);
    }
    // cout << "mesonAmplitude = " << mesonAmplitude << endl;
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;

    // special setting for PCM-EMC
    // if (icalo == 2 || icalo == 14){
    mesonAmplitudeMin = mesonAmplitude*80./100.;
    mesonAmplitudeMax = mesonAmplitude*110./100.;
    if (icalo == 1){
      fitRangeMax = 0.2;
    }
    // cout << "mesonAmplitudeMin = " << mesonAmplitudeMin << endl;
    // cout << "mesonAmplitudeMax = " << mesonAmplitudeMax << endl;

    TF1* fFitReco    = new TF1("fGaussExp","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                               fitRangeMin, fitRangeMax);
    Double_t fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();

    fFitReco->SetParameter(0, mesonAmplitude);
    fFitReco->SetParameter(1, fMesonMassExpect);
    fFitReco->SetParameter(2, 0.01);
    fFitReco->SetParameter(3, 0.012);
    fFitReco->SetParameter(5, 0.0);
    fFitReco->SetParLimits(5, 0.000001,0.000002);
    fFitReco->SetParLimits(4, 0.0,1000);

    fFitReco->SetParLimits(0, mesonAmplitudeMin, mesonAmplitudeMax);
    // if (icalo == 4 || icalo == 12 || icalo == 15){
      fFitReco->SetParLimits(1, fMesonMassExpect*0.75, fMesonMassExpect*1.2);
    // } else {
        // fFitReco->SetParLimits(1, fMesonMassExpect*0.9, fMesonMassExpect*1.1);
    // }
    fFitReco->SetParLimits(2, 0.001, 0.1);
    fFitReco->SetParLimits(3, 0.001, 0.09);
    // if (icalo == 5){
    //   fFitReco->SetParLimits(1, fMesonMassExpect*0.8, fMesonMassExpect*1.8);
    //   fFitReco->SetParLimits(2, 0.001, 0.01);
    // }

    histo->Fit(fFitReco,"QRME0");
    histo->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(kRed+1);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    // if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
    //     cout << "Parameter for exponential+Gaussian "<< endl;
    //     cout << gMinuit->fCstatu.Data() << endl;
    //     cout << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<<endl;
    //     cout << "chi2/ndf:" << fFitReco->GetChisquare()/fFitReco->GetNDF() << endl;
    // } else {
    //     cout << "Exp+Gaussian fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
    // }
    return fFitReco;
}

TF1* FitExpPlusGaussianPol2(TH1F* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t icalo, Double_t ptcenter , Int_t iDataMC){

    Double_t mesonAmplitude = 0;
    for(Int_t i = histo->FindBin(0.09); i < histo->FindBin(0.2) ; i++ ){
      if(histo->GetBinContent(i) > mesonAmplitude) mesonAmplitude = histo->GetBinContent(i);
    }
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;

    mesonAmplitudeMin = mesonAmplitude*70./100.;
    mesonAmplitudeMax = mesonAmplitude*110./100.;
    // if(icalo==3){
    //   mesonAmplitudeMin = mesonAmplitude*50./100.;
    // }
    if (icalo == 2){
      fitRangeMin = 0.07;
      fitRangeMax = 0.17;
    } else if (icalo == 1){
      // fitRangeMax = 0.19;
    }

    TF1* fFitReco    = new TF1("GaussExpPol2","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",fitRangeMin, fitRangeMax);
    // TF1* fFitReco    = new TF1("GaussExpPol2","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x+[6]*x*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x+[6]*x*x)",fitRangeMin, fitRangeMax);
    Double_t fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    // amplitude
    fFitReco->SetParameter(0, mesonAmplitude);
    fFitReco->SetParLimits(0, mesonAmplitudeMin, mesonAmplitudeMax);

    // peak position
    fFitReco->SetParameter(1, fMesonMassExpect);
    fFitReco->SetParLimits(1, fMesonMassExpect*0.75, fMesonMassExpect*1.1);

    // peak width
    fFitReco->SetParameter(2, 0.03);
    fFitReco->SetParLimits(2, 0.01, 0.5);
    // if(icalo==0){
    //   fFitReco->SetParLimits(2, 0.001, 0.03);
    // }

    // lambda tail
    fFitReco->SetParameter(3, 0.012);
    // fFitReco->SetParLimits(3, 0.00001, 0.09);

    // constant offset
    fFitReco->SetParLimits(4, 0.0,1000);

    // linear BG
    fFitReco->SetParameter(5, 1.0);

    // pol2 BG
    // fFitReco->SetParameter(6, 0.1);
    // fFitReco->SetParLimits(6, 0.0, 0.0);



    histo->Fit(fFitReco,"QRME0");
    histo->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(kRed+1);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    return fFitReco;
}
