#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 

void eoverpstudies(
    TString inputFileName     = "central_allculsters.root",
    TString suffix            = "pdf",
    Bool_t doFitting          = kTRUE,
    TString collisionsSys     = "PythiaMB"
){  
  TString clusterizerName = "MA clusters";
  const int maxcalo = 3;
  TString caloName[maxcalo+1] = {"FEMC", "EEMC", "BECAL", "Hybrid"};
  bool caloactive[maxcalo] = {false};
  Double_t mass_pi0PDG = 0.1349770;

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


  
  TString outputDir                 = Form("plotsEoverP");
  if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFile                      = new TFile(inputFileName.Data());
  if(inputFile->IsZombie()){
    cout << "inputFile not found!" << endl;
  }

  //************************** Read data **************************************************
  const int nEne = 15;
  const static Double_t partE[]   = {   0.7, 1.0, 1.5, 2, 3, 4, 6, 10,15,20,25,  30,35,  40, 50};
  double energymax = 40.;
  // const int nEne = 26;
  // const static Double_t partE[]   = {   0.7, 0.9, 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 4.5,  5,
  //                                       5.5, 6, 7, 8, 10,
  //                                      12.5,15,17.5,20,22.5,25,  30,35,  40, 50};
  // double energymax = 40.;
  int padebin = 0;
  int padebinelec = 0;
  TH1F*  h_TMstudies_EoverP_trackSpec_pions[maxcalo] = {NULL}; // [calorimeter_enum][algorithm_enum]
  TH1F*  h_TMstudies_EoverP_trackSpec_electrons[maxcalo] = {NULL}; // [calorimeter_enum][algorithm_enum]
  TH2F*  h_TMstudies_EoverP_pions_EoverP_E[maxcalo] = {NULL}; // [calorimeter_enum][algorithm_enum]
  TH2F*  h_TMstudies_EoverP_electrons_EoverP_E[maxcalo] = {NULL}; // [calorimeter_enum][algorithm_enum]

  TH1F* histEoP_proj_pions[maxcalo][50] = {{NULL}};
  TH1F* histEoP_proj_electrons[maxcalo][50] = {{NULL}};
  TH1F* histIntCounts_pions_EoP[maxcalo] = {NULL};
  TH1F* histIntCounts_electrons_EoP[maxcalo] = {NULL};
  TH1F* histIntCounts_pions_all[maxcalo] = {NULL};
  TH1F* histIntCounts_electrons_all[maxcalo] = {NULL};
  TH1F* histPionRejection[maxcalo] = {NULL};


  for(int icalo=0;icalo<maxcalo;icalo++){
    cout << "processing " << caloName[icalo].Data() << endl;
    TString NameEP_hist_EpE = Form("%s/h_TMstudies_EoverP_pions_EoverP_E_%s",caloName[icalo].Data(), caloName[icalo].Data() );
    // TString NameE_Minv = Form("%s/h_Etrue_Minv_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    h_TMstudies_EoverP_pions_EoverP_E[icalo]    = (TH2F*)inputFile->Get(NameEP_hist_EpE.Data());
    if(!h_TMstudies_EoverP_pions_EoverP_E[icalo]){
      cout << "\t" <<NameEP_hist_EpE.Data() << " not found" << endl;
      continue;
    } else {
      caloactive[icalo]=true;
    };
    h_TMstudies_EoverP_electrons_EoverP_E[icalo]    = (TH2F*)inputFile->Get(Form("%s/h_TMstudies_EoverP_electrons_EoverP_E_%s",caloName[icalo].Data(), caloName[icalo].Data() ));
    h_TMstudies_EoverP_trackSpec_pions[icalo]    = (TH1F*)inputFile->Get(Form("%s/h_TMstudies_EoverP_trackSpec_pions_%s",caloName[icalo].Data(), caloName[icalo].Data() ));
    h_TMstudies_EoverP_trackSpec_electrons[icalo]    = (TH1F*)inputFile->Get(Form("%s/h_TMstudies_EoverP_trackSpec_electrons_%s",caloName[icalo].Data(), caloName[icalo].Data() ));


//     cout << "\tinputs loaded" << endl;

    if (h_TMstudies_EoverP_pions_EoverP_E[icalo]->GetEntries() > 100 || h_TMstudies_EoverP_electrons_EoverP_E[icalo]->GetEntries()>100){
      float latextopedge = 0.90;
      float latexsideedge = 0.92;
      float  textSizeLabelsRel      = 55./1200;
      TCanvas* Minvpi02D = new TCanvas(Form("Minvpi02D%s",caloName[icalo].Data()),"",0,0,1500,800);
      DrawGammaCanvasSettings( Minvpi02D, 0.095, 0.01, 0.01, 0.105);
      Minvpi02D->Divide(2,1);
      Minvpi02D->cd(1);
      h_TMstudies_EoverP_pions_EoverP_E[icalo]->GetXaxis()->SetRangeUser(0.,2);
      // h_TMstudies_EoverP_pions_EoverP_E[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      SetStyleHistoTH2ForGraphs(h_TMstudies_EoverP_pions_EoverP_E[icalo],"E^{MC}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      h_TMstudies_EoverP_pions_EoverP_E[icalo]->Draw("colz");
      drawLatexAdd("true #pi^{-} tracks",0.17,latextopedge,textSizeLabelsRel,false,false,false);      
      Minvpi02D->cd(2);
      h_TMstudies_EoverP_electrons_EoverP_E[icalo]->GetXaxis()->SetRangeUser(0.,2);
      // h_TMstudies_EoverP_electrons_EoverP_E[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      SetStyleHistoTH2ForGraphs(h_TMstudies_EoverP_electrons_EoverP_E[icalo],"E^{reco}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
      h_TMstudies_EoverP_electrons_EoverP_E[icalo]->Draw("colz");
      
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
      Minvpi02D->Print(Form("%s/EoverP_vs_E_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));
    }
    padebin = 0;
    padebinelec = 0;

  // TH1F* histEoP_proj_pions[maxcalo][50] = {{NULL}};
  // TH1F* histEoP_proj_electrons[maxcalo][50] = {{NULL}};

    histIntCounts_pions_EoP[icalo] = new TH1F(Form("histIntCounts_pions_EoP_%s",caloName[icalo].Data()),"",nEne-1,partE);
    histIntCounts_electrons_EoP[icalo] = new TH1F(Form("histIntCounts_electrons_EoP_%s",caloName[icalo].Data()),"",nEne-1,partE);
    histIntCounts_pions_all[icalo] = new TH1F(Form("histIntCounts_pions_all_%s",caloName[icalo].Data()),"",nEne-1,partE);
    histIntCounts_electrons_all[icalo] = new TH1F(Form("histIntCounts_electrons_all_%s",caloName[icalo].Data()),"",nEne-1,partE);
    histPionRejection[icalo] = new TH1F(Form("histPionRejection_%s",caloName[icalo].Data()),"",nEne-1,partE);
    // histMeanEVsMinv[icalo] = new TH1F(Form("histMeanEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
    // histSigmaEVsMinv[icalo] = new TH1F(Form("histSigmaEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
    for(int ebin=0;ebin<nEne-1;ebin++){
      histEoP_proj_pions[icalo][padebin]    = (TH1F*)h_TMstudies_EoverP_pions_EoverP_E[icalo]->ProjectionX(Form("histEoP_proj_pions_%s_%f",caloName[icalo].Data(), partE[ebin]), 
                                            h_TMstudies_EoverP_pions_EoverP_E[icalo]->GetYaxis()->FindBin(partE[ebin]), h_TMstudies_EoverP_pions_EoverP_E[icalo]->GetYaxis()->FindBin(partE[ebin+1]),"e");
      histEoP_proj_electrons[icalo][padebinelec]    = (TH1F*)h_TMstudies_EoverP_electrons_EoverP_E[icalo]->ProjectionX(Form("histEoP_proj_electrons_%s_%f",caloName[icalo].Data(), partE[ebin]), 
                                            h_TMstudies_EoverP_electrons_EoverP_E[icalo]->GetYaxis()->FindBin(partE[ebin]), h_TMstudies_EoverP_electrons_EoverP_E[icalo]->GetYaxis()->FindBin(partE[ebin+1]),"e");

      if(histEoP_proj_pions[icalo][padebin]->GetEntries() > 10){
        // TH1F* histIntCounts_pions_EoP[maxcalo] = {NULL};
        // TH1F* histIntCounts_electrons_EoP[maxcalo] = {NULL};
        // TH1F* histIntCounts_pions_all[maxcalo] = {NULL};
        // TH1F* histIntCounts_electrons_all[maxcalo] = {NULL};

        Double_t intpions  = histEoP_proj_pions[icalo][padebin]->Integral(histEoP_proj_pions[icalo][padebin]->GetXaxis()->FindBin(0.7), histEoP_proj_pions[icalo][padebin]->GetXaxis()->FindBin(1.2));
        histIntCounts_pions_EoP[icalo]->SetBinContent(histIntCounts_pions_EoP[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpions);
        
        Double_t intpionsall  = h_TMstudies_EoverP_trackSpec_pions[icalo]->Integral(h_TMstudies_EoverP_trackSpec_pions[icalo]->GetXaxis()->FindBin(partE[ebin]), h_TMstudies_EoverP_trackSpec_pions[icalo]->GetXaxis()->FindBin(partE[ebin]+1));
        histIntCounts_pions_all[icalo]->SetBinContent(histIntCounts_pions_all[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpionsall);
        // if(icalo==2){
        //   histEoP_proj_pions[icalo][padebin]->Rebin(6);
        //   // histEoP_proj_pions_allcls[icalo][padebin]->Rebin(6);
        // } else {
        //   histEoP_proj_pions[icalo][padebin]->Rebin(4);
        //   // histEoP_proj_pions_allcls[icalo][padebin]->Rebin(4);
        // }
        // mean = meanErr = sigma = sigmaErr = 0;

        // if (doFitting){   

        //   fithistEtrueMinv[icalo][padebin] = FitExpPlusGaussian (histEoP_proj_pions[icalo][padebin], minmassfit, maxmassfit, icalo, (partE[ebin]+partE[ebin+1])/2, 0);
        //   // fithistEtrueMinv[icalo][padebin]    = new TF1(Form("fit_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]), "crystalball", 0.05, 0.15);
        //   // fithistEtrueMinv[icalo][padebin]->SetParameters(0.55*histEoP_proj_pions[icalo][padebin]->GetMaximum(),histEoP_proj_pions[icalo][padebin]->GetMean(),2.*histEoP_proj_pions[icalo][padebin]->GetRMS(),2.*histEoP_proj_pions[icalo][padebin]->GetRMS(),2.5*histEoP_proj_pions[icalo][padebin]->GetRMS());
        //   histEoP_proj_pions[icalo][padebin]->Fit(fithistEtrueMinv[icalo][padebin],"L0RMEQ","",minmassfit, maxmassfit);

        //   // fithistEtrueMinvFINAL[icalo][padebin] = histEoP_proj_pions[icalo][padebin]->GetFunction(Form("fit_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]));
        //   fithistEtrueMinvFINAL[icalo][padebin] = fithistEtrueMinv[icalo][padebin];

        //   mean     = fithistEtrueMinvFINAL[icalo][padebin]->GetParameter(1);
        //   meanErr  = fithistEtrueMinvFINAL[icalo][padebin]->GetParError(1);
        //   sigma    = fithistEtrueMinvFINAL[icalo][padebin]->GetParameter(2);
        //   sigmaErr = fithistEtrueMinvFINAL[icalo][padebin]->GetParError(2);
        //   allmean[icalo].push_back(mean);
        //   allsigma[icalo].push_back(sigma);
        // }


        // Double_t reso     = sigma/mean;
        // Double_t resoErr  = TMath::Sqrt( TMath::Power(sigmaErr/sigma,2) + TMath::Power(meanErr/mean,2) );

        // histResoEVsMinv[icalo]->SetBinContent(histResoEVsMinv[icalo]->FindBin(partE[ebin]),  100*reso);
        // histResoEVsMinv[icalo]->SetBinError(histResoEVsMinv[icalo]->FindBin(partE[ebin]),  100*resoErr);  
        // histMeanEVsMinv[icalo]->SetBinContent(histMeanEVsMinv[icalo]->FindBin(partE[ebin]), mean);
        // histMeanEVsMinv[icalo]->SetBinError(histMeanEVsMinv[icalo]->FindBin(partE[ebin]), meanErr);

        // histSigmaEVsMinv[icalo]->SetBinContent(histSigmaEVsMinv[icalo]->FindBin(partE[ebin]), sigma/sigmat);
        // histSigmaEVsMinv[icalo]->SetBinError(histSigmaEVsMinv[icalo]->FindBin(partE[ebin]), sigmaErr/sigmat);
        padebin++;
      }
      if(histEoP_proj_electrons[icalo][padebinelec]->GetEntries() > 10){
        Double_t intelectrons  = histEoP_proj_electrons[icalo][padebinelec]->Integral(histEoP_proj_electrons[icalo][padebinelec]->GetXaxis()->FindBin(0.7), histEoP_proj_electrons[icalo][padebinelec]->GetXaxis()->FindBin(1.2));
        histIntCounts_electrons_EoP[icalo]->SetBinContent(histIntCounts_electrons_EoP[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectrons);

        Double_t intelectronsall  = h_TMstudies_EoverP_trackSpec_electrons[icalo]->Integral(h_TMstudies_EoverP_trackSpec_electrons[icalo]->GetXaxis()->FindBin(partE[ebin]), h_TMstudies_EoverP_trackSpec_electrons[icalo]->GetXaxis()->FindBin(partE[ebin]+1));
        histIntCounts_electrons_all[icalo]->SetBinContent(histIntCounts_electrons_all[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectronsall);

        padebinelec++;
      }
      cout << (partE[ebin]+partE[ebin+1])/2 << "\tpions: " << histIntCounts_pions_EoP[icalo]->GetBinContent(histIntCounts_pions_EoP[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << "\telectrons: " << histIntCounts_electrons_EoP[icalo]->GetBinContent(histIntCounts_electrons_EoP[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << endl;
      cout << "\tallpions: " << histIntCounts_pions_all[icalo]->GetBinContent(histIntCounts_pions_all[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << "\tallelectrons: " << histIntCounts_electrons_all[icalo]->GetBinContent(histIntCounts_electrons_all[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << endl;
      double rejection = (histIntCounts_pions_all[icalo]->GetBinContent(histIntCounts_pions_all[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) / histIntCounts_electrons_all[icalo]->GetBinContent(histIntCounts_electrons_all[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2))) / (histIntCounts_pions_EoP[icalo]->GetBinContent(histIntCounts_pions_EoP[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) / histIntCounts_electrons_EoP[icalo]->GetBinContent(histIntCounts_electrons_EoP[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)));
      if(isnan(rejection)) rejection = -1;
      if(isinf(rejection)) rejection = -1;
      cout << "\t\trejection: " << rejection << endl;
      histPionRejection[icalo]->SetBinContent(histPionRejection[icalo]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),rejection);
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
        Tl.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", caloName[icalo].Data()));
        Tl.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
        Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

        legend1->Draw();
      } else{
        mass2D->cd(z+1);
        if(!histEoP_proj_pions[icalo][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
        //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
        SetStyleHistoTH1ForGraphs( histEoP_proj_pions[icalo][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
        DrawGammaSetMarker(histEoP_proj_pions[icalo][padnumber],34,1, kBlack, kBlack);
        histEoP_proj_pions[icalo][padnumber]->GetYaxis()->SetNdivisions(5);
        histEoP_proj_pions[icalo][padnumber]->GetXaxis()->SetNdivisions(5);
        histEoP_proj_pions[icalo][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
        histEoP_proj_pions[icalo][padnumber]->Draw("ep");
        // DrawGammaSetMarker(histEoP_proj_pions_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
        // histEoP_proj_pions_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
        // histEoP_proj_pions_allcls[icalo][padnumber]->Draw("same,hist");
        // if(fithistEtrueMinvFINAL[icalo][padnumber]){
        //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
        //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
        //   fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
        //   DrawGammaLines(allmean[icalo][padnumber], allmean[icalo][padnumber], 0, 0.5*fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(),2, 2);
        //   DrawGammaLines(allmean[icalo][padnumber] + allsigma[icalo][padnumber], allmean[icalo][padnumber]+ allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
        //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
        //   if(z==0){
        //     legend1->AddEntry(histEoP_proj_pions[icalo][padnumber], "Simulation", "p");
        //     legend1->AddEntry(fithistEtrueMinvFINAL[icalo][padnumber], "Fit", "l");
        //   }
      drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
      }
      padnumber++;
    }

    mass2D->Print(Form("%s/EoP_bin_pions%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));

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
        Tl.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", caloName[icalo].Data()));
        Tl.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
        Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

        legend1->Draw();
      } else{
        mass2D->cd(z+1);
        if(!histEoP_proj_electrons[icalo][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
        //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
        SetStyleHistoTH1ForGraphs( histEoP_proj_electrons[icalo][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
        DrawGammaSetMarker(histEoP_proj_electrons[icalo][padnumber],34,1, kBlack, kBlack);
        histEoP_proj_electrons[icalo][padnumber]->GetYaxis()->SetNdivisions(5);
        histEoP_proj_electrons[icalo][padnumber]->GetXaxis()->SetNdivisions(5);
        histEoP_proj_electrons[icalo][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
        histEoP_proj_electrons[icalo][padnumber]->Draw("ep");
        // DrawGammaSetMarker(histEoP_proj_electrons_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
        // histEoP_proj_electrons_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
        // histEoP_proj_electrons_allcls[icalo][padnumber]->Draw("same,hist");
        // if(fithistEtrueMinvFINAL[icalo][padnumber]){
        //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
        //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
        //   fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
          DrawGammaLines(0.7, 0.7, 0, 0.5*histEoP_proj_electrons[icalo][padnumber]->GetMaximum(),1, kGray+1,3);
          DrawGammaLines(1.2, 1.2, 0, 0.5*histEoP_proj_electrons[icalo][padnumber]->GetMaximum(),1, kGray+1,3);
        //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
        //   if(z==0){
        //     legend1->AddEntry(histEoP_proj_pions[icalo][padnumber], "Simulation", "p");
        //     legend1->AddEntry(fithistEtrueMinvFINAL[icalo][padnumber], "Fit", "l");
        //   }
      drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
      }
      padnumber++;
    }

    mass2D->Print(Form("%s/EoP_bin_electrons%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));

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
        Tl.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", caloName[icalo].Data()));
        Tl.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
        Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

        legend1->Draw();
      } else{
        mass2D->cd(z+1);
        if(!histEoP_proj_electrons[icalo][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
        //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
        SetStyleHistoTH1ForGraphs( histEoP_proj_electrons[icalo][padnumber], "#it{p}_{track}", "counts", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,0.6,1.1);
        DrawGammaSetMarker(histEoP_proj_electrons[icalo][padnumber],34,1, kBlack, kBlack);
        histEoP_proj_electrons[icalo][padnumber]->GetYaxis()->SetNdivisions(5);
        histEoP_proj_electrons[icalo][padnumber]->GetXaxis()->SetNdivisions(5);
        histEoP_proj_electrons[icalo][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
        histEoP_proj_electrons[icalo][padnumber]->Draw("ep");
        DrawGammaSetMarker(histEoP_proj_pions[icalo][padnumber],28,1, kRed+2, kRed+2);
        // histEoP_proj_pions[icalo][padnumber]->GetYaxis()->SetNdivisions(5);
        // histEoP_proj_pions[icalo][padnumber]->GetXaxis()->SetNdivisions(5);
        // histEoP_proj_pions[icalo][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
        histEoP_proj_pions[icalo][padnumber]->Draw("ep,same");
        // DrawGammaSetMarker(histEoP_proj_electrons_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
        // histEoP_proj_electrons_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
        // histEoP_proj_electrons_allcls[icalo][padnumber]->Draw("same,hist");
        // if(fithistEtrueMinvFINAL[icalo][padnumber]){
        //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
        //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
        //   fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
          DrawGammaLines(0.7, 0.7, 0, 0.5*histEoP_proj_electrons[icalo][padnumber]->GetMaximum(),1, kGray+1,3);
          DrawGammaLines(1.2, 1.2, 0, 0.5*histEoP_proj_electrons[icalo][padnumber]->GetMaximum(),1, kGray+1,3);
        //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
          if(z==0){
            legend1->AddEntry(histEoP_proj_electrons[icalo][padnumber], "Electron E/p", "p");
            legend1->AddEntry(histEoP_proj_pions[icalo][padnumber], "Pion E/p", "p");
          }
      drawLatexAdd(Form("%1.1f < #it{p}_{track} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
      }
      padnumber++;
    }

    mass2D->Print(Form("%s/EoP_bin_electrons_and_pions_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));
  }



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
  histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.51, 49, 1000, 1, 999 );
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{track} (GeV)", "#pi^{#pm}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510,42,42);//(#times #epsilon_{pur})
  histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
  histo2DAccEff->GetXaxis()->SetNoExponent();
  histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAccEff->DrawCopy();
  TGraphAsymmErrors* graphPionRejection[maxcalo]={NULL};
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
      if(histPionRejection[icalo]){
        graphPionRejection[icalo] = new TGraphAsymmErrors(histPionRejection[icalo]);
        while(graphPionRejection[icalo]->GetY()[graphPionRejection[icalo]->GetN()-1] == -1) graphPionRejection[icalo]->RemovePoint(graphPionRejection[icalo]->GetN()-1);
        DrawGammaSetMarkerTGraphAsym(graphPionRejection[icalo], markerStylePID[icalo], 1.2*markerSizePID[icalo], colorPID[icalo] , colorPID[icalo]);
        graphPionRejection[icalo]->Draw("samelpe");
      }
  }

  TLegend* legendPi0Merge           = GetAndSetLegend2(0.15, 0.92-0.1, 0.3, 0.92-0.1-(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  for (Int_t icalo = maxcalo-1; icalo > -1; icalo--){
      if(histPionRejection[icalo]){
          legendPi0Merge->AddEntry(graphPionRejection[icalo],caloName[icalo].Data(),"p");
      }
  }
  legendPi0Merge->Draw();
  float toplegval = 0.90;
  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  // DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


  canvasMergingFraction->Update();
  canvasMergingFraction->Print(Form("%s/EoP_pionrejection.%s",outputDir.Data(),suffix.Data()));
 


}

//       if(!allcalofitted){
//         // while(ebin < energymax){
//         int padebin = 0;

//         histResoEVsMinv_allcalo = new TH1F(Form("histResoEVsMinv_%s","hybrid"),"",nEne-1,partE);
//         histMeanEVsMinv_allcalo = new TH1F(Form("histMeanEVsMinv_%s","hybrid"),"",nEne-1,partE);
//         histSigmaEVsMinv_allcalo = new TH1F(Form("histSigmaEVsMinv_%s","hybrid"),"",nEne-1,partE);
//         for(int ebin=0;ebin<nEne-1;ebin++){
//           histEtrueMinvbin_allcalo[padebin]    = (TH1F*)histEtrueMinv_allcalo->ProjectionY(Form("projection_Minv_%s_%f","hybrid", partE[ebin]), 
//                                                 histEtrueMinv_allcalo->GetXaxis()->FindBin(partE[ebin]), histEtrueMinv_allcalo->GetXaxis()->FindBin(partE[ebin+1]),"e");

//           if(histEtrueMinvbin_allcalo[padebin]->GetEntries() > 10){
//             histEtrueMinvbin_allcalo[padebin]->Rebin(4);
//             mean = meanErr = sigma = sigmaErr = 0;

//             if (doFitting){

//               // fithistEtrueMinv_allcalo[padebin]    = new TF1(Form("fit_Minv_%s_%f","hybrid", partE[ebin]), "crystalball", 0.05, 0.15);
//               // fithistEtrueMinv_allcalo[padebin]->SetParameters(0.55*histEtrueMinvbin_allcalo[padebin]->GetMaximum(),histEtrueMinvbin_allcalo[padebin]->GetMean(),2.*histEtrueMinvbin_allcalo[padebin]->GetRMS(),2.*histEtrueMinvbin_allcalo[padebin]->GetRMS(),2.5*histEtrueMinvbin_allcalo[padebin]->GetRMS());
//               fithistEtrueMinv_allcalo[padebin] = FitExpPlusGaussian (histEtrueMinvbin_allcalo[padebin], minmassfit, maxmassfit, icalo, (partE[ebin]+partE[ebin+1])/2, 0);
//               histEtrueMinvbin_allcalo[padebin]->Fit(fithistEtrueMinv_allcalo[padebin],"L0RMEQ","",minmassfit, maxmassfit);
              
//               fithistEtrueMinvFINAL_allcalo[padebin] = fithistEtrueMinv_allcalo[padebin];
//               // fithistEtrueMinvFINAL_allcalo[padebin] = histEtrueMinvbin_allcalo[padebin]->GetFunction(Form("fit_Minv_%s_%f","hybrid", partE[ebin]));
                  
//               mean     = fithistEtrueMinvFINAL_allcalo[padebin]->GetParameter(1);
//               meanErr  = fithistEtrueMinvFINAL_allcalo[padebin]->GetParError(1);
//               sigma    = fithistEtrueMinvFINAL_allcalo[padebin]->GetParameter(2);
//               sigmaErr = fithistEtrueMinvFINAL_allcalo[padebin]->GetParError(2);
//               allmean_allcalo.push_back(mean);
//               allsigma_allcalo.push_back(sigma);
//             }


//             Double_t reso     = sigma/mean;
//             Double_t resoErr  = TMath::Sqrt( TMath::Power(sigmaErr/sigma,2) + TMath::Power(meanErr/mean,2) );
              
//             histResoEVsMinv_allcalo->SetBinContent(histResoEVsMinv_allcalo->FindBin(partE[ebin]),  100*reso);
//             histResoEVsMinv_allcalo->SetBinError(histResoEVsMinv_allcalo->FindBin(partE[ebin]),  100*resoErr);  
//             histMeanEVsMinv_allcalo->SetBinContent(histMeanEVsMinv_allcalo->FindBin(partE[ebin]), mean);
//             histMeanEVsMinv_allcalo->SetBinError(histMeanEVsMinv_allcalo->FindBin(partE[ebin]), meanErr);

//             histSigmaEVsMinv_allcalo->SetBinContent(histSigmaEVsMinv_allcalo->FindBin(partE[ebin]), sigma/sigmat);
//             histSigmaEVsMinv_allcalo->SetBinError(histSigmaEVsMinv_allcalo->FindBin(partE[ebin]), sigmaErr/sigmat);
//             padebin++;
//           }
//           // ebin = ebin + 2.; 
//         }


//         mass2D = new TCanvas("mass2D","",0,0,1000,800);
//         gStyle->SetTitleX(0.5);
//         gStyle->SetTitleAlign(23);
//         gStyle->SetTitleBorderSize(0);
//         gStyle->SetTitleFontSize(0.07); 
//         DrawGammaCanvasSettings( mass2D, 0.095, 0.01, 0.01, 0.10);
//         mass2D->Divide(5,6);
//         int padnumber = 0;
//         legend1 = new TLegend(0., 0.01, 0.4, 0.25);
//         for(int z=0; z< padebin +1 ; z++ ){
//           if(z==4){
//             mass2D->cd(z+1);
//             TLatex Tl;
//             Tl.SetTextSize(0.1);
//             Tl.DrawLatex(0.03,0.9,perfLabel);
//             Tl.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", "Hybrid"));
//             Tl.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
//             Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

//             legend1->Draw();
//           } else{
//             mass2D->cd(z+1);
//             //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
//             SetStyleHistoTH1ForGraphs( histEtrueMinvbin_allcalo[padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
//             DrawGammaSetMarker(histEtrueMinvbin_allcalo[padnumber],34,1, kBlack, kBlack);
//             histEtrueMinvbin_allcalo[padnumber]->GetYaxis()->SetNdivisions(5);
//             histEtrueMinvbin_allcalo[padnumber]->GetXaxis()->SetNdivisions(5);
//             histEtrueMinvbin_allcalo[padnumber]->GetXaxis()->SetRangeUser(0,0.3);
//             histEtrueMinvbin_allcalo[padnumber]->Draw("ep");
//             // DrawGammaSetMarker(histEtrueMinvbin_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
//             // histEtrueMinvbin_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
//             // histEtrueMinvbin_allcls[icalo][padnumber]->Draw("same,hist");
//             fithistEtrueMinvFINAL_allcalo[padnumber]->SetLineColor(kBlue);
//             fithistEtrueMinvFINAL_allcalo[padnumber]->SetLineWidth(1);
//             fithistEtrueMinvFINAL_allcalo[padnumber]->Draw("same");
//             drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
//             DrawGammaLines(allmean_allcalo[padnumber], allmean_allcalo[padnumber], 0, 0.5*fithistEtrueMinvFINAL_allcalo[padnumber]->GetMaximum(),2, 2);
//             DrawGammaLines(allmean_allcalo[padnumber] + allsigma_allcalo[padnumber], allmean_allcalo[padnumber]+ allsigma_allcalo[padnumber], 0, fithistEtrueMinvFINAL_allcalo[padnumber]->GetMaximum(), 1, 12);
//             DrawGammaLines(allmean_allcalo[padnumber] - allsigma_allcalo[padnumber], allmean_allcalo[padnumber]- allsigma_allcalo[padnumber], 0, fithistEtrueMinvFINAL_allcalo[padnumber]->GetMaximum(), 1, 12);
//             if(z==0){
//               legend1->AddEntry(histEtrueMinvbin_allcalo[padnumber], "Simulation", "p");
//               legend1->AddEntry(fithistEtrueMinvFINAL_allcalo[padnumber], "Fit", "l");
//             }
//             padnumber++;
//           }
//         }
//         mass2D->Print(Form("%s/Mass_%s_%s.%s", outputDir.Data(), "hybrid", clusterizerName.Data(), suffix.Data()));

//         allcalofitted = true;
//       }

//       //============= Fit Results======================================
//       TCanvas* fitresults = new TCanvas("fitresults","",0,0,1000,500);
//       DrawGammaCanvasSettings( fitresults, 0.15, 0.15, 0.15, 0.15);
    
//       fitresults->Divide(2,1);
//       fitresults->cd(1);

//       histMeanEVsMinv[icalo]->GetYaxis()->SetRangeUser(0.08,0.16);
//       SetStyleHistoTH1ForGraphs( histMeanEVsMinv[icalo], "(M_{#gamma #gamma}", "", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.,1.5);
//       DrawGammaSetMarker(histMeanEVsMinv[icalo],34,1, kBlack, kBlack);   
//       histMeanEVsMinv[icalo]->GetYaxis()->SetRangeUser(0,1.3*histMeanEVsMinv[icalo]->GetMaximum()); 
//       histMeanEVsMinv[icalo]->Draw("e,p");
//       DrawGammaLines(0, energymax, mass_pi0PDG, mass_pi0PDG,2, 13);

//       TLatex Tl3;
//       Tl3.SetTextSize(0.035);
//       Tl3.DrawLatex(2,0.18,perfLabel);
//       Tl3.DrawLatex(2,0.17,Form("Calorimeter: %s ", caloName[icalo].Data()));
//       Tl3.DrawLatex(2,0.16,Form("Clusterization type: %s ", clusterizerName.Data()));
//       Tl3.DrawLatex(2,0.14,"#pi ^{0} #rightarrow #gamma #gamma");
//       Tl3.Draw();
//       fitresults->cd(2);

//       histSigmaEVsMinv[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
//       SetStyleHistoTH1ForGraphs(histSigmaEVsMinv[icalo], "E^{MC}", "Width (M_{#gamma #gamma})", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.3);  
//       DrawGammaSetMarker(histSigmaEVsMinv[icalo],34,1, kBlack, kBlack);   
//       histSigmaEVsMinv[icalo]->GetYaxis()->SetRangeUser(0,0.015); 
//       histSigmaEVsMinv[icalo]->Draw("e,p");
      
      
//       fitresults->Print(Form("%s/Mean_Width_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));
//     }


//     TString Name_mclus = Form("%s/h_Pi0Spec_mergedclus_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
//     h_Pi0Spec_mergedclus[icalo]    = (TH1F*)inputFile->Get(Name_mclus.Data());
//     TString Name_sepclus = Form("%s/h_Pi0Spec_separatecls_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
//     h_Pi0Spec_separatecls[icalo]    = (TH1F*)inputFile->Get(Name_sepclus.Data());
//     TString Name_totclus = Form("%s/h_Pi0Spec_total_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
//     h_Pi0Spec_total[icalo]    = (TH1F*)inputFile->Get(Name_totclus.Data());
//   }

  
//   //============= End Fit Results======================================

//   float minPtPi0 = 0.6;
//   float maxPtPi0woMerged = 49;

//   Double_t mesonMassExpectPi0                     = TDatabasePDG::Instance()->GetParticle(111)->Mass();
//   Double_t mesonMassExpectEta                     = TDatabasePDG::Instance()->GetParticle(221)->Mass();
//   Double_t arrayBoundariesX1_4[2];
//   Double_t arrayBoundariesY1_4[3];
//   Double_t relativeMarginsX[3];
//   Double_t relativeMarginsY[3];
//   textSizeLabelsPixel             = 50;
//   ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

//   TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
//   DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

//   TPad* padWidthPi0               = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
//   DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
//   padWidthPi0->Draw();

//   TPad* padMassPi0                = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
//   DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
//   padMassPi0->Draw();

//   TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
//   DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
//   padMassLegend1->SetFillStyle(0);
//   padMassLegend1->Draw();

//   padWidthPi0->cd();
//   padWidthPi0->SetLogx();

//   Double_t margin                 = relativeMarginsX[0]*2.7*1350;
//   Double_t textsizeLabelsWidth    = 0;
//   Double_t textsizeFacWidth       = 0;
//   if (padWidthPi0->XtoPixel(padWidthPi0->GetX2()) < padWidthPi0->YtoPixel(padWidthPi0->GetY1())){
//       textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
//       textsizeFacWidth            = (Double_t)1./padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
//   } else {
//       textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->YtoPixel(padWidthPi0->GetY1());
//       textsizeFacWidth            = (Double_t)1./padWidthPi0->YtoPixel(padWidthPi0->GetY1());
//   }

//   TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, minPtPi0,maxPtPi0woMerged ,1000., -30, 40);
//   SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{E} (GeV)", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
//                             0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin),512,505,42,42);//#it{p}_{T} (GeV/#it{c})
//   histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(0.1,19.5);//24.5);
//   histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
//   histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
//   histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
//   histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
//   histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
//   histo2DAllPi0FWHM->DrawCopy();

//   Color_t colorDet[maxcalo+1]                          = {kBlack, kRed+2, kBlue+2, kGreen+2};//, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2,kOrange+2,kSpring+2};
//   Style_t markerStyleDet[maxcalo+1]   = {20, 47, 28, 34};//, 28, 30, 42, 46, 24, 25, 27, 28, 30};
//   Size_t markerSizeDet[maxcalo+1]     = {2.5, 3, 3 , 3};//, 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9,1.5, 1.8 };
//   TGraphAsymmErrors* graphPi0FWHMMeV[maxcalo];
//   for (Int_t icalo = 0; icalo < maxcalo; icalo++){
//     if(caloactive[icalo]){
//       if(histSigmaEVsMinv[icalo]){
//         histSigmaEVsMinv[icalo]->Scale(1000);
//         graphPi0FWHMMeV[icalo] = new TGraphAsymmErrors(histSigmaEVsMinv[icalo]);
//         if(graphPi0FWHMMeV[icalo]){
//           if(icalo==1){
//             while(graphPi0FWHMMeV[icalo]->GetX()[graphPi0FWHMMeV[icalo]->GetN()-1] > 30) graphPi0FWHMMeV[icalo]->RemovePoint(graphPi0FWHMMeV[icalo]->GetN()-1);
//           } else if(icalo==2){
//             while(graphPi0FWHMMeV[icalo]->GetX()[graphPi0FWHMMeV[icalo]->GetN()-1] > 12) graphPi0FWHMMeV[icalo]->RemovePoint(graphPi0FWHMMeV[icalo]->GetN()-1);
//           }
//             DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[icalo], markerStyleDet[icalo], markerSizeDet[icalo], colorDet[icalo] , colorDet[icalo]);
//             // DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[0] , colorDet[0]);
//             graphPi0FWHMMeV[icalo]->Draw("p,same,e");
//             // DrawGammaSetMarkerTGraphAsym(graphPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
//             // graphPi0TrueFWHMMeV[i]->Draw("p,same,e");
//         }
//       }
//     }
//   }
//   TGraphAsymmErrors* graphPi0FWHMMeV_allcalo;
//   if(histSigmaEVsMinv_allcalo){
//     histSigmaEVsMinv_allcalo->Scale(1000);
//     graphPi0FWHMMeV_allcalo = new TGraphAsymmErrors(histSigmaEVsMinv_allcalo);
//     if(graphPi0FWHMMeV_allcalo){
//       while(graphPi0FWHMMeV_allcalo->GetX()[graphPi0FWHMMeV_allcalo->GetN()-1] > 30) graphPi0FWHMMeV_allcalo->RemovePoint(graphPi0FWHMMeV_allcalo->GetN()-1);
//         DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV_allcalo, markerStyleDet[maxcalo], markerSizeDet[maxcalo], colorDet[maxcalo] , colorDet[maxcalo]);
//         // DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[0] , colorDet[0]);
//         graphPi0FWHMMeV_allcalo->Draw("p,same,e");
//         // DrawGammaSetMarkerTGraphAsym(graphPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
//         // graphPi0TrueFWHMMeV[i]->Draw("p,same,e");
//     }
//   }

//   // TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
//   // SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
//   // labelLegendAMass->SetTextFont(43);
//   // labelLegendAMass->Draw();

//   TLatex *labelMassPerf       = new TLatex(0.13,0.87,perfLabel.Data());//"ALICE performance");
//   SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
//   labelMassPerf->SetTextFont(43);
//   labelMassPerf->Draw();
//   TLatex *labelMassEnergy     = new TLatex(0.13,0.78,Form("%s clusters", clusterizerName.Data()));
//   SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4);
//   labelMassEnergy->SetTextFont(43);
//   labelMassEnergy->Draw();
//   TLatex *labelMassPi0        = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
//   SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4);
//   labelMassPi0->SetTextFont(43);
//   labelMassPi0->Draw();

//   padMassPi0->cd();
//   padMassPi0->SetLogx();

//   Double_t textsizeLabelsMass         = 0;
//   Double_t textsizeFacMass            = 0;
//   if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
//       textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
//       textsizeFacMass                 = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
//   } else {
//       textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
//       textsizeFacMass                 = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
//   }

//   TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minPtPi0,maxPtPi0woMerged, 1000., 80.1, 167.9);//, 100.1, 160.9);//125.1, 155.9);
//   SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{E}_{#pi^{0}} (GeV)", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
//                             textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin),512,505,42,42);
//   histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
//   histo2DAllPi0Mass->GetYaxis()->SetRangeUser(80.1, 149.9);//(131.1, 147.9);//125.1, 155.9);
//   histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
//   histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
//   histo2DAllPi0Mass->GetXaxis()->SetNoExponent();
//   histo2DAllPi0Mass->DrawCopy();
//   DrawGammaLines(minPtPi0,maxPtPi0woMerged , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,2, kGray+1,7);
//   TGraphAsymmErrors* graphPi0Mass[maxcalo];
//   TLegend* legendMassWidth   = GetAndSetLegend2(0.12, 0.2, 0.32, 0.2+4*textsizeLabelsMass, textSizeLabelsPixel);
//   for (Int_t icalo = 0; icalo < maxcalo; icalo++){
//     if(caloactive[icalo] && histMeanEVsMinv[icalo]){
//       histMeanEVsMinv[icalo]->Scale(1000);
//       graphPi0Mass[icalo] = new TGraphAsymmErrors(histMeanEVsMinv[icalo]);
//       if(icalo==1){
//           while(graphPi0Mass[icalo]->GetX()[graphPi0Mass[icalo]->GetN()-1] > 30) graphPi0Mass[icalo]->RemovePoint(graphPi0Mass[icalo]->GetN()-1);
//         } else if(icalo==2){
//           while(graphPi0Mass[icalo]->GetX()[graphPi0Mass[icalo]->GetN()-1] > 12) graphPi0Mass[icalo]->RemovePoint(graphPi0Mass[icalo]->GetN()-1);
//         }
//       if(graphPi0Mass[icalo]){
//         DrawGammaSetMarkerTGraphAsym(graphPi0Mass[icalo], markerStyleDet[icalo], markerSizeDet[icalo], colorDet[icalo] , colorDet[icalo]);
//         graphPi0Mass[icalo]->Draw("p,same,e");
//         // DrawGammaSetMarkerTGraphAsym(graphPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
//         // graphPi0TrueMass[i]->Draw("p,same,e");
//         legendMassWidth->AddEntry(graphPi0Mass[icalo],caloName[icalo].Data(),"p");
//       }
//     }
//   }
//   TGraphAsymmErrors* graphPi0Mass_allcalo;
//   if(histMeanEVsMinv_allcalo){
//     histMeanEVsMinv_allcalo->Scale(1000);
//     graphPi0Mass_allcalo = new TGraphAsymmErrors(histMeanEVsMinv_allcalo);
//     if(graphPi0Mass_allcalo){
//       while(graphPi0Mass_allcalo->GetX()[graphPi0Mass_allcalo->GetN()-1] > 30) graphPi0Mass_allcalo->RemovePoint(graphPi0Mass_allcalo->GetN()-1);
//       DrawGammaSetMarkerTGraphAsym(graphPi0Mass_allcalo, markerStyleDet[maxcalo], markerSizeDet[maxcalo], colorDet[maxcalo] , colorDet[maxcalo]);
//       graphPi0Mass_allcalo->Draw("p,same,e");
//       // DrawGammaSetMarkerTGraphAsym(graphPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
//       // graphPi0TrueMass[i]->Draw("p,same,e");
//       legendMassWidth->AddEntry(graphPi0Mass_allcalo,"Hybrid","p");
//     }
//   }
//   legendMassWidth->Draw();


//   // TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
//   // SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
//   // labelLegendBMass->SetTextFont(43);
//   // labelLegendBMass->Draw();

//   // legendM02Data->SetMargin(0.05/(0.9-0.8));
//   // //********************************** Defintion of the Legend **************************************************
//   // Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
//   // Double_t  rowsLegendMass2[14]= {0.84,0.66,0.50,0.33,0.16,0.00,0.16,0.16,0.16,0.16,0.16,0.16,0.16};
//   // //******************* Offsets ***********************
//   // Double_t offsetMarkerXMass2         = 0.1;
//   // Double_t offsetMarkerYMass2         = 0.1;
//   // //****************** Scale factors ******************
//   // Double_t scaleMarkerMass2           = 1.2;

//   // padMassLegend1->cd();
//   // //****************** first Column **************************************************
//   // TLatex *textMassPCM[10];
//   // Int_t legendRunningIndex = 0;
//   // for (Int_t icalo = 0; icalo < maxcalo; icalo++){
//   //     if(graphPi0Mass[icalo]){
//   //         textMassPCM[icalo]                  = new TLatex(columnsLegendMass2[0],rowsLegendMass2[legendRunningIndex+1],caloName[icalo].Data());
//   //         SetStyleTLatex( textMassPCM[icalo], textSizeLabelsPixel,4);
//   //         textMassPCM[icalo]->SetTextFont(43);
//   //         textMassPCM[icalo]->Draw();
//   //         legendRunningIndex++;
//   //     }
//   // }
//   // //****************** second Column *************************************************
//   // TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
//   // SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
//   // textMassData->SetTextFont(43);
//   // textMassData->Draw();
//   // TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
//   // SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
//   // textMassMC->SetTextFont(43);
//   // textMassMC->Draw();

//   // TMarker* markerPCMPi0Mass[10];
//   // TMarker* markerPCMPi0MassMC[10];
//   // legendRunningIndex = 0;
//   // for (Int_t i = 0; i < 13; i++){
//   //     if(graphPi0Mass[i] && graphPi0TrueMass[i]){
//   //         markerPCMPi0Mass[i]             = CreateMarkerFromGraph(graphPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
//   //         markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2);
//   //         markerPCMPi0MassMC[i]           = CreateMarkerFromGraph(graphPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
//   //         markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2);
//   //         legendRunningIndex++;
//   //     }
//   // }

//   canvasMassWidthPi0->Update();
//   canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));




//     textSizeLabelsPixel                 = 100*3.5/5;
//     TCanvas* canvasInvMassSamplePlot    = new TCanvas("canvasInvMassSamplePlot","",0,0,1500,1500);  // gives the page size
//     DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.105, 0.01, 0.039, 0.093);
//     Double_t textSizeLabelPPMass             = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;

//     Style_t markerStyleInvMassSGBG      = 0;
//     Size_t markerSizeInvMassSGBG        = 0;
//     Color_t markerColorInvMassSGBG      = kBlack;
//     Style_t markerStyleInvMassMBG       = 24;
//     Size_t markerSizeInvMassMBG         = 3;
//     Color_t markerColorInvMassMBG       = kGray+2;
//     Color_t markerColorInvMassMBG1      = kGray+3;
//     Color_t markerColorInvMassMBG2      = kGray+1;
//     Style_t markerStyleInvMassBG        = 20;
//     Size_t markerSizeInvMassBG          = 3;
//     Color_t markerColorInvMassBG        = kBlack;
//     Style_t markerStyleInvMassSG        = 20;
//     Size_t markerSizeInvMassSG          = 4;
//     Color_t markerColorInvMassSG        = kRed+2;
//     Color_t fitColorInvMassSG           = kAzure+2;

//     Double_t marginInvMass          = 0.1*1500;
//     Double_t textsizeLabelsInvMass  = 0;
//     Double_t textsizeFacInvMass     = 0;
//     if (canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) < canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1())){
//         textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
//         textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
//     } else {
//         textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
//         textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
//     }
//     cout << textsizeLabelsInvMass << endl;

//     TH2F * histo2DPi0InvMassDummy;
//     histo2DPi0InvMassDummy             = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.05,0.21,21000,-1000,200000);
//     // histo2DPi0InvMassDummy             = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.05,0.249,21000,-1000,200000);
//     SetStyleHistoTH2ForGraphs(histo2DPi0InvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
//                             0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass),505,510,42,42);
//   int exampleBins[maxcalo+1] = {16,7,6,11};
//   bool plotinvmass = false;
//   if(plotinvmass){
//     for (Int_t icalo = 0; icalo < maxcalo+1; icalo++){
//       if(histEtrueMinvbin[icalo][exampleBins[icalo]]){
//         Int_t iptbinplot = exampleBins[icalo];
//         TH1D* histoSignalPlusBGPi0                        = nullptr;
//         if(icalo<maxcalo){
//           histoSignalPlusBGPi0                        = (TH1D*)histEtrueMinvbin[icalo][iptbinplot]->Clone(Form("Mapping_GG_InvMass_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
//         } else {
//           histoSignalPlusBGPi0                        = (TH1D*)histEtrueMinvbin_allcalo[iptbinplot]->Clone(Form("Mapping_GG_InvMass_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
//         }
//         TH1D* histoSignalPi0                              = nullptr;
//         TH1D* histoCombBGPi0                              = nullptr;
//         TF1* fFitGausExp                                  = nullptr;
//         if(icalo<maxcalo){
//           fFitGausExp  = (TF1*)fithistEtrueMinvFINAL[icalo][iptbinplot]->Clone(Form("Signal_InvMassFit_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
//         } else {
//           fFitGausExp = (TF1*)fithistEtrueMinvFINAL_allcalo[iptbinplot]->Clone(Form("Signal_InvMassFit_in_Pt_Bin%d_%s",iptbinplot,caloName[icalo].Data()));
//         }
//         fFitGausExp->SetNpx(10000);

//         canvasInvMassSamplePlot->cd();
//         // histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
//         histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(1.2*histoSignalPlusBGPi0->GetMinimum(),1.2*histoSignalPlusBGPi0->GetMaximum());
//         histo2DPi0InvMassDummy->DrawCopy();

//         TLatex *labelInvMassPtRangel = new TLatex(0.945,0.9,Form("#pi^{0}: %2.1f < #it{E}_{#pi^{0}} < %2.1f GeV/#it{c}",partE[iptbinplot],partE[iptbinplot+1]));
//         // TLatex *labelInvMassPtRangel = new TLatex(0.945,0.9,Form("#pi^{0}: %2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",partE[iptbinplot],partE[iptbinplot+1]));
//         if(histoSignalPlusBGPi0){
//           // DrawGammaSetMarker(histoSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
//           DrawGammaSetMarker(histoSignalPlusBGPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSGBG, markerColorInvMassSGBG);
//           histoSignalPlusBGPi0->SetLineWidth(1);
//           // histoSignalPlusBGPi0->Draw("hist,e,same");
//           histoSignalPlusBGPi0->Draw("same");
//         }
//         // if(histoCombBGPi0){
//         //   DrawGammaSetMarker(histoCombBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
//         //   histoCombBGPi0->SetLineWidth(markerSizeInvMassMBG);
//         //   histoCombBGPi0->Draw("same");
//         // }
//         // if(histoSignalPi0){
//         //   DrawGammaSetMarker(histoSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
//         //   histoSignalPi0->SetStats(kFALSE);
//         //   histoSignalPi0->Draw("same");
//         // }
//         fFitGausExp->SetLineWidth(1);
//         fFitGausExp->SetRange(0,0.255);
//         fFitGausExp->SetLineColor(fitColorInvMassSG);
//         fFitGausExp->Draw("same");

//         TLatex *labelInvMassALICE      = new TLatex(0.135,0.9,perfLabel.Data());
//         SetStyleTLatex( labelInvMassALICE, 0.85*textSizeLabelsPixel,4);
//         labelInvMassALICE->SetTextFont(43);
//         labelInvMassALICE->Draw();

//         // TLatex *labelInvMassInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textSizeLabelPPMass,collisionSystempp8TeV.Data());
//         TLatex *labelInvMassInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.9*textSizeLabelPPMass,collisionSystem.Data());
//         SetStyleTLatex( labelInvMassInvMassEnergy, 0.85*textSizeLabelsPixel,4);
//         labelInvMassInvMassEnergy->SetTextFont(43);
//         labelInvMassInvMassEnergy->Draw();

//         TLatex *labelInvMassInvMassTriggerl      = new TLatex(0.135,0.9-0.9*2*0.9*textSizeLabelPPMass,caloName[icalo].Data());
//         SetStyleTLatex( labelInvMassInvMassTriggerl, 0.85*textSizeLabelsPixel,4);
//         labelInvMassInvMassTriggerl->SetTextFont(43);
//         labelInvMassInvMassTriggerl->Draw();

//         // TLatex *labelInvMassInvMassReco  = new TLatex(0.135,0.9-0.9*3*0.8*textSizeLabelPPMass,"EMC");
//         // SetStyleTLatex( labelInvMassInvMassReco, 0.85*textSizeLabelsPixel,4);
//         // labelInvMassInvMassReco->SetTextFont(43);
//         // labelInvMassInvMassReco->Draw();

//         SetStyleTLatex( labelInvMassPtRangel, 0.85*textSizeLabelsPixel,4);
//         labelInvMassPtRangel->SetTextAlign(31);
//         labelInvMassPtRangel->SetTextFont(43);
//         labelInvMassPtRangel->Draw();

//         DrawGammaLines(fFitGausExp->GetParameter(1)-(0.06),fFitGausExp->GetParameter(1)-(0.06),1.2*histoSignalPlusBGPi0->GetMinimum(),0.2*histoSignalPlusBGPi0->GetMaximum(),2,kGray+2,7);
//         DrawGammaLines(fFitGausExp->GetParameter(1)+(0.080),fFitGausExp->GetParameter(1)+(0.080),1.2*histoSignalPlusBGPi0->GetMinimum(),0.2*histoSignalPlusBGPi0->GetMaximum(),2,kGray+2,7);
//         // TLegend* legendInvMassl  =GetAndSetLegend2(0.62, 0.87-5*0.8*0.75*textSizeLabelPPMass, 0.85, 0.87, 0.85*textSizeLabelsPixel);
//         int legendreduction = 0;
//         if(!histoCombBGPi0)legendreduction+=2;
//         if(!histoSignalPi0)legendreduction++;
//         TLegend* legendInvMassl  = GetAndSetLegend2(0.62, 0.87-(5-legendreduction)*0.9*0.75*textSizeLabelPPMass, 0.85, 0.87, 0.85*textSizeLabelsPixel);
//         legendInvMassl->SetMargin(0.25);
//         legendInvMassl->AddEntry(histoSignalPlusBGPi0,"Raw real events","p");
//         // legendInvMassl->AddEntry(histoSignalPlusBGPi0,"Raw real events","l");
//         if(histoCombBGPi0){
//           legendInvMassl->AddEntry(histoCombBGPi0,"Mixed event +","p");
//           legendInvMassl->AddEntry((TObject*)0,"remain. BG","");
//         } else {
//           // legendInvMassl->AddEntry((TObject*)0,"Mixed event +","");
//           // legendInvMassl->AddEntry((TObject*)0,"remain. BG","");
//         }
//         if(histoSignalPi0){
//           legendInvMassl->AddEntry(histoSignalPi0,"BG subtracted","p");
//         } else {
//           // legendInvMassl->AddEntry((TObject*)0,"BG subtracted","");
//         }
//         legendInvMassl->AddEntry(fFitGausExp, "Fit","l");
//         legendInvMassl->Draw();
//         canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBin%d_%s.%s",outputDir.Data(),iptbinplot,caloName[icalo].Data(), suffix.Data()));
//       }
//     }
//   }


//     // h_Pi0Spec_mergedclus[icalo]    = (TH2F*)inputFile->Get(Name_mclus.Data());
//     // h_Pi0Spec_separatecls[icalo]    = (TH2F*)inputFile->Get(Name_sepclus.Data());
//     // h_Pi0Spec_total[icalo]    = (TH2F*)inputFile->Get(Name_totclus.Data());
//   TH1F* h_Pi0_mergedfraction[maxcalo] = {NULL};
//   TH1F* h_Pi0_separatefraction[maxcalo] = {NULL};

//   // **********************************************************************************************************************
//   // ******************************** Acceptance * Efficiency for pi0 single measurement 8TeV **************************
//   // **********************************************************************************************************************
//   textSizeLabelsPixel             = 55;
//   textSizeLabelsRel      = 55./1200;
//   cout << textSizeLabelsRel << endl;

//   TCanvas* canvasMergingFraction       = new TCanvas("canvasMergingFraction", "", 200, 10, 1200, 1100);  // gives the page size
//   DrawGammaCanvasSettings( canvasMergingFraction,  0.1, 0.01, 0.015, 0.095);
//   // canvasMergingFraction->SetLogy(1);
//   // canvasMergingFraction->SetLogx(1);

//   TH2F * histo2DAccEff;
//   histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, minPtPi0, maxPtPi0woMerged, 1000, 0, 1.09 );
//   SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{E}_{#pi^{0}} (GeV)", "#it{f}_{merg.}",
//                           0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510,42,42);//(#times #epsilon_{pur})
//   histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
//   histo2DAccEff->GetXaxis()->SetNoExponent();
//   histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histo2DAccEff->DrawCopy();
//   for (Int_t icalo = 0; icalo < maxcalo; icalo++){
//       if(h_Pi0Spec_mergedclus[icalo]){
//         h_Pi0Spec_mergedclus[icalo]->Sumw2();
//         h_Pi0Spec_mergedclus[icalo]->Rebin(4);
//         h_Pi0Spec_total[icalo]->Sumw2();
//         h_Pi0Spec_total[icalo]->Rebin(4);
//         h_Pi0_mergedfraction[icalo] = (TH1F*)h_Pi0Spec_mergedclus[icalo]->Clone(Form("mergedfrachisto_%d",icalo));
//         // h_Pi0_mergedfraction[icalo]->Divide(h_Pi0Spec_total[icalo]);
//         h_Pi0_mergedfraction[icalo]->Divide(h_Pi0_mergedfraction[icalo],h_Pi0Spec_total[icalo], 1, 1,"B" );
//           DrawGammaSetMarker(h_Pi0_mergedfraction[icalo], markerStylePID[icalo], 1.2*markerSizePID[icalo], colorPID[icalo] , colorPID[icalo]);
//           h_Pi0_mergedfraction[icalo]->Draw("same,p,e");
//       }
//   }

//   TLegend* legendPi0Merge           = GetAndSetLegend2(0.45, 0.13, 0.95, 0.13+(1*textSizeLabelsRel),0.9*textSizeLabelsPixel,3);
//   for (Int_t icalo = maxcalo-1; icalo > -1; icalo--){
//       if(h_Pi0Spec_mergedclus[icalo]){
//           legendPi0Merge->AddEntry(h_Pi0_mergedfraction[icalo],caloName[icalo].Data(),"p");
//       }
//   }
//   legendPi0Merge->Draw();
//   float lowlegval = 0.20;
//   drawLatexAdd(perfLabel.Data(),0.95,lowlegval+0.1,textSizeLabelsRel,false,false,true);
//   drawLatexAdd(collisionSystem.Data(),0.95,lowlegval+0.05,textSizeLabelsRel,false,false,true);
//   drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.95,lowlegval,textSizeLabelsRel,false,false,true);

//   DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


//   canvasMergingFraction->Update();
//   canvasMergingFraction->Print(Form("%s/Pi0_MergingFraction.%s",outputDir.Data(),suffix.Data()));
 
// }


// TF1* FitExpPlusGaussian(TH1F* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t icalo, Double_t ptcenter , Int_t iDataMC){

//     Double_t mesonAmplitude = 0;
//     for(Int_t i = histo->FindBin(0.07); i < histo->FindBin(0.2) ; i++ ){
//       if(histo->GetBinContent(i) > mesonAmplitude) mesonAmplitude = histo->GetBinContent(i);
//     }
//     // cout << "mesonAmplitude = " << mesonAmplitude << endl;
//     Double_t mesonAmplitudeMin;
//     Double_t mesonAmplitudeMax;

//     // special setting for PCM-EMC
//     // if (icalo == 2 || icalo == 14){
//     mesonAmplitudeMin = mesonAmplitude*80./100.;
//     mesonAmplitudeMax = mesonAmplitude*110./100.;
//     if (icalo == 2){
//       fitRangeMax = 0.2;
//     }
//     // cout << "mesonAmplitudeMin = " << mesonAmplitudeMin << endl;
//     // cout << "mesonAmplitudeMax = " << mesonAmplitudeMax << endl;

//     TF1* fFitReco    = new TF1("fGaussExp","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
//                                fitRangeMin, fitRangeMax);
//     Double_t fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();

//     fFitReco->SetParameter(0, mesonAmplitude);
//     fFitReco->SetParameter(1, fMesonMassExpect);
//     fFitReco->SetParameter(2, 0.01);
//     fFitReco->SetParameter(3, 0.012);
//     fFitReco->SetParameter(5, 0.0);
//     fFitReco->SetParLimits(5, 0.000001,0.000002);
//     fFitReco->SetParLimits(4, 0.0,1000);

//     fFitReco->SetParLimits(0, mesonAmplitudeMin, mesonAmplitudeMax);
//     // if (icalo == 4 || icalo == 12 || icalo == 15){
//         fFitReco->SetParLimits(1, fMesonMassExpect*0.75, fMesonMassExpect*1.2);
//     // } else {
//         // fFitReco->SetParLimits(1, fMesonMassExpect*0.9, fMesonMassExpect*1.1);
//     // }
//     fFitReco->SetParLimits(2, 0.001, 0.1);
//     fFitReco->SetParLimits(3, 0.001, 0.09);
//     // if (icalo == 5){
//     //   fFitReco->SetParLimits(1, fMesonMassExpect*0.8, fMesonMassExpect*1.8);
//     //   fFitReco->SetParLimits(2, 0.001, 0.01);
//     // }

//     histo->Fit(fFitReco,"QRME0");
//     histo->Fit(fFitReco,"QRME0");

//     fFitReco->SetLineColor(kRed+1);
//     fFitReco->SetLineWidth(1);
//     fFitReco->SetLineStyle(1);

//     // if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
//     //     cout << "Parameter for exponential+Gaussian "<< endl;
//     //     cout << gMinuit->fCstatu.Data() << endl;
//     //     cout << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<<endl;
//     //     cout << "chi2/ndf:" << fFitReco->GetChisquare()/fFitReco->GetNDF() << endl;
//     // } else {
//     //     cout << "Exp+Gaussian fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
//     // }
//     return fFitReco;
// }