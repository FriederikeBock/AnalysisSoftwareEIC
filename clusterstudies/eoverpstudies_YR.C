#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 

void eoverpstudies_YR(
    TString inputFileNamePion     = "central_allculsters.root",
    TString inputFileNameElectron     = "central_allculsters.root",
    TString suffix            = "pdf",
    Bool_t doFitting          = kTRUE,
    TString collisionsSys     = "PythiaMB"
){  
  TString clusterizerName = "MA clusters";
  const int maxcalo = 3;
  TString caloName[maxcalo+1] = {"FEMC", "EEMC", "BECAL", "Hybrid"};
  bool caloactive[maxcalo] = {false};
  Double_t mass_pi0PDG = 0.1349770;
  double eopcutvalue[maxcalo] = {0.8, 0.75, 0.77};
  double eopcutvaluemax[maxcalo] = {2.0, 2.0, 2.0};
  double etarangemin[maxcalo] = {1.2, -4, -1.6};
  double etarangemax[maxcalo] = {4.0, -1.4, 1.2};
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
  double energymax = 40.;
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
  TH3F*  h_EoverP_pions_EoverP_E_Eta[maxmatch][maxcalo] = {{NULL}}; // [calorimeter_enum][algorithm_enum]
  TH3F*  h_EoverP_electrons_EoverP_E_Eta[maxmatch][maxcalo] = {{NULL}}; // [calorimeter_enum][algorithm_enum]

  TH1F* histEoP_proj_pions[maxmatch][maxcalo][nEta][50]     = {{{{NULL}}}};
  TH1F* histEoP_proj_electrons[maxmatch][maxcalo][nEta][50] = {{{{NULL}}}};
  TH1F* histIntCounts_pions_EoP[maxmatch][maxcalo][nEta]    = {{{NULL}}};
  TH1F* histIntCounts_electrons_EoP[maxmatch][maxcalo][nEta]= {{{NULL}}};
  TH1F* histIntCounts_pions_all[maxmatch][maxcalo][nEta]    = {{{NULL}}};
  TH1F* histIntCounts_electrons_all[maxmatch][maxcalo][nEta]= {{{NULL}}};
  TH1F* histPionRejection[maxmatch][maxcalo][nEta]          = {{{NULL}}};


  // TString matchbasis = "trackbased";
  TString matchbasis[maxmatch] = {"trackbased", "projectionbased"};
  for(int imatch=0;imatch<maxmatch;imatch++){
    for(int icalo=0;icalo<maxcalo;icalo++){
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

      h_EoverP_pions_EoverP_E_Eta[imatch][icalo]    = (TH3F*)inputFilePion->Get(Form("%s/h_EoverP_pions_EoverP_E_Eta_%s_%s",caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data() ));
      h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]    = (TH3F*)inputFileElec->Get(Form("%s/h_EoverP_electrons_EoverP_E_Eta_%s_%s",caloName[icalo].Data(), matchbasis[imatch].Data(), caloName[icalo].Data() ));




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
        Minvpi02D->Print(Form("%s/EoverP_vs_E_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));
      }


      for(int ieta=0;ieta<nEta;ieta++){
        cout << partEta[ieta] << " < eta < " << partEta[ieta+1] << endl;
        padebin = 0;
        padebinelec = 0;
        histIntCounts_pions_EoP[imatch][icalo][ieta] = new TH1F(Form("histIntCounts_pions_EoP_%s_%f",caloName[icalo].Data(),partEta[ieta]),"",nEne-1,partE);
        histIntCounts_electrons_EoP[imatch][icalo][ieta] = new TH1F(Form("histIntCounts_electrons_EoP_%s_%f",caloName[icalo].Data(),partEta[ieta]),"",nEne-1,partE);
        histIntCounts_pions_all[imatch][icalo][ieta] = new TH1F(Form("histIntCounts_pions_all_%s_%f",caloName[icalo].Data(),partEta[ieta]),"",nEne-1,partE);
        histIntCounts_electrons_all[imatch][icalo][ieta] = new TH1F(Form("histIntCounts_electrons_all_%s_%f",caloName[icalo].Data(),partEta[ieta]),"",nEne-1,partE);
        histPionRejection[imatch][icalo][ieta] = new TH1F(Form("histPionRejection_%s_%f",caloName[icalo].Data(),partEta[ieta]),"",nEne-1,partE);


        for(int ebin=0;ebin<nEne-1;ebin++){
          h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->GetYaxis()->SetRange(h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->GetYaxis()->FindBin(partE[ebin]),h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->GetYaxis()->FindBin(partE[ebin+1]));
          h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->GetZaxis()->SetRange(h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->GetZaxis()->FindBin(partEta[ieta]),h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->GetZaxis()->FindBin(partEta[ieta+1]));
          h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->GetYaxis()->SetRange(h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->GetYaxis()->FindBin(partE[ebin]),h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->GetYaxis()->FindBin(partE[ebin+1]));
          h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->GetZaxis()->SetRange(h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->GetZaxis()->FindBin(partEta[ieta]),h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->GetZaxis()->FindBin(partEta[ieta+1]));

          histEoP_proj_pions[imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->Project3D(Form("%1.1f_%1.1f_x",partEta[ieta],partE[ebin]));
          histEoP_proj_electrons[imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->Project3D(Form("%1.1f_%1.1f_x",partEta[ieta],partE[ebin]));
          // histEoP_proj_pions[imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_pions_EoverP_E_Eta[imatch][icalo]->Project3D(Form("histEoP_proj_pions_%s_%f_%f",caloName[icalo].Data(),partEta[ieta], partE[ebin]));
          // histEoP_proj_electrons[imatch][icalo][ieta][padebin]    = (TH1F*)h_EoverP_electrons_EoverP_E_Eta[imatch][icalo]->Project3D(Form("histEoP_proj_pions_%s_%f_%f",caloName[icalo].Data(),partEta[ieta], partE[ebin]));


          Double_t intpions  = histEoP_proj_pions[imatch][icalo][ieta][padebin]->Integral(histEoP_proj_pions[imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvalue[icalo]), histEoP_proj_pions[imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          histIntCounts_pions_EoP[imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_EoP[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpions);

          Double_t intpionsall  = histEoP_proj_pions[imatch][icalo][ieta][padebin]->Integral(histEoP_proj_pions[imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(0.0), histEoP_proj_pions[imatch][icalo][ieta][padebin]->GetXaxis()->FindBin(1.9));
          histIntCounts_pions_all[imatch][icalo][ieta]->SetBinContent(histIntCounts_pions_all[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intpionsall);
          padebin++;


          Double_t intelectrons  = histEoP_proj_electrons[imatch][icalo][ieta][padebinelec]->Integral(histEoP_proj_electrons[imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(eopcutvalue[icalo]), histEoP_proj_electrons[imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(eopcutvaluemax[icalo]));
          histIntCounts_electrons_EoP[imatch][icalo][ieta]->SetBinContent(histIntCounts_electrons_EoP[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectrons);

          Double_t intelectronsall  = histEoP_proj_electrons[imatch][icalo][ieta][padebinelec]->Integral(histEoP_proj_electrons[imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(0.0), histEoP_proj_electrons[imatch][icalo][ieta][padebinelec]->GetXaxis()->FindBin(1.9));
          histIntCounts_electrons_all[imatch][icalo][ieta]->SetBinContent(histIntCounts_electrons_all[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),intelectronsall);

          padebinelec++;

          cout << (partE[ebin]+partE[ebin+1])/2 << "\tpions: " << histIntCounts_pions_EoP[imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_EoP[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << "\telectrons: " << histIntCounts_electrons_EoP[imatch][icalo][ieta]->GetBinContent(histIntCounts_electrons_EoP[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << endl;
          cout << "\tallpions: " << histIntCounts_pions_all[imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_all[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << "\tallelectrons: " << histIntCounts_electrons_all[imatch][icalo][ieta]->GetBinContent(histIntCounts_electrons_all[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) << endl;
          double rejection = histIntCounts_pions_all[imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_all[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2)) / histIntCounts_pions_EoP[imatch][icalo][ieta]->GetBinContent(histIntCounts_pions_EoP[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2));
          if(isnan(rejection)) rejection = -1;
          if(isinf(rejection)) rejection = -1;
          cout << "\t\trejection: " << rejection << endl;
          histPionRejection[imatch][icalo][ieta]->SetBinContent(histPionRejection[imatch][icalo][ieta]->GetXaxis()->FindBin((partE[ebin]+partE[ebin+1])/2),rejection);
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
          Tl.DrawLatex(0.03,0.45,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.3,Form("%1.2f < #eta < %1.1f", partEta[ieta], partEta[ieta+1]));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          if(!histEoP_proj_pions[imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEoP_proj_pions[imatch][icalo][ieta][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
          DrawGammaSetMarker(histEoP_proj_pions[imatch][icalo][ieta][padnumber],34,1, kBlack, kBlack);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->Draw("ep");
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
          //     legend1->AddEntry(histEoP_proj_pions[imatch][icalo][ieta][padnumber], "Simulation", "p");
          //     legend1->AddEntry(fithistEtrueMinvFINAL[icalo][padnumber], "Fit", "l");
          //   }
        drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        }
        padnumber++;
      }
      if(partEta[ieta] > etarangemin[icalo] && partEta[ieta+1] < etarangemax[icalo])
        mass2D->Print(Form("%s/EoP_bin_pions%s_%s_%f.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), partEta[ieta], suffix.Data()));

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
          Tl.DrawLatex(0.03,0.45,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.3,Form("%1.2f < #eta < %1.1f", partEta[ieta], partEta[ieta+1]));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          if(!histEoP_proj_electrons[imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEoP_proj_electrons[imatch][icalo][ieta][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
          DrawGammaSetMarker(histEoP_proj_electrons[imatch][icalo][ieta][padnumber],34,1, kBlack, kBlack);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->Draw("ep");
          // DrawGammaSetMarker(histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
          // histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
          // histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber]->Draw("same,hist");
          // if(fithistEtrueMinvFINAL[icalo][padnumber]){
          //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
          //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
          //   fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
            DrawGammaLines(eopcutvalue[icalo], eopcutvalue[icalo], 0, 0.5*histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            DrawGammaLines(eopcutvaluemax[icalo], eopcutvaluemax[icalo], 0, 0.5*histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
          //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
          //   if(z==0){
          //     legend1->AddEntry(histEoP_proj_pions[imatch][icalo][ieta][padnumber], "Simulation", "p");
          //     legend1->AddEntry(fithistEtrueMinvFINAL[icalo][padnumber], "Fit", "l");
          //   }
        drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        }
        padnumber++;
      }

      if(partEta[ieta] > etarangemin[icalo] && partEta[ieta+1] < etarangemax[icalo])
      mass2D->Print(Form("%s/EoP_bin_electrons%s_%s_%f.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), partEta[ieta], suffix.Data()));

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
          Tl.DrawLatex(0.03,0.45,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.3,Form("%1.2f < #eta < %1.1f", partEta[ieta], partEta[ieta+1]));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          if(!histEoP_proj_pions[imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEoP_proj_pions[imatch][icalo][ieta][padnumber], "#it{p}_{track}", "counts", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,0.6,1.1);
          DrawGammaSetMarker(histEoP_proj_pions[imatch][icalo][ieta][padnumber],34,1, kRed+2, kRed+2);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->Draw("ep");
            DrawGammaLines(eopcutvalue[icalo], eopcutvalue[icalo], 0, 0.5*histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            DrawGammaLines(eopcutvaluemax[icalo], eopcutvaluemax[icalo], 0, 0.5*histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
          //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
            if(z==0){
              legend1->AddEntry(histEoP_proj_pions[imatch][icalo][ieta][padnumber], "Pion E/p", "p");
            }
        drawLatexAdd(Form("%1.1f < #it{p}_{track} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        }
        padnumber++;
      }

      if(partEta[ieta] > etarangemin[icalo] && partEta[ieta+1] < etarangemax[icalo])
      mass2D->Print(Form("%s/EoP_bin_pions_%s_%s_%f.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), partEta[ieta], suffix.Data()));

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
          Tl.DrawLatex(0.03,0.45,"#pi ^{0} #rightarrow #gamma #gamma");
          Tl.DrawLatex(0.03,0.3,Form("%1.2f < #eta < %1.1f", partEta[ieta], partEta[ieta+1]));

          legend1->Draw();
        } else{
          mass2D->cd(z+1);
          if(!histEoP_proj_electrons[imatch][icalo][ieta][padnumber]){ cout << "lol" << endl;padnumber++; continue;}
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEoP_proj_electrons[imatch][icalo][ieta][padnumber], "#it{p}_{track}", "counts", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,0.6,1.1);
          DrawGammaSetMarker(histEoP_proj_electrons[imatch][icalo][ieta][padnumber],34,1, kBlack, kBlack);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->Draw("ep");
          DrawGammaSetMarker(histEoP_proj_pions[imatch][icalo][ieta][padnumber],28,1, kRed+2, kRed+2);
          // histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetYaxis()->SetNdivisions(5);
          // histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetXaxis()->SetNdivisions(5);
          // histEoP_proj_pions[imatch][icalo][ieta][padnumber]->GetXaxis()->SetRangeUser(0,2.0);
          histEoP_proj_pions[imatch][icalo][ieta][padnumber]->Draw("ep,same");
          // DrawGammaSetMarker(histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
          // histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
          // histEoP_proj_electrons[imatch][icalo][ieta]_allcls[icalo][padnumber]->Draw("same,hist");
          // if(fithistEtrueMinvFINAL[icalo][padnumber]){
          //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
          //   fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
          //   fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
            DrawGammaLines(eopcutvalue[icalo], eopcutvalue[icalo], 0, 0.5*histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
            DrawGammaLines(eopcutvaluemax[icalo], eopcutvaluemax[icalo], 0, 0.5*histEoP_proj_electrons[imatch][icalo][ieta][padnumber]->GetMaximum(),1, kGray+1,3);
          //   DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
            if(z==0){
              legend1->AddEntry(histEoP_proj_electrons[imatch][icalo][ieta][padnumber], "Electron E/p", "p");
              legend1->AddEntry(histEoP_proj_pions[imatch][icalo][ieta][padnumber], "Pion E/p", "p");
            }
        drawLatexAdd(Form("%1.1f < #it{p}_{track} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
        }
        padnumber++;
      }

      if(partEta[ieta] > etarangemin[icalo] && partEta[ieta+1] < etarangemax[icalo])
      mass2D->Print(Form("%s/EoP_bin_electrons_and_pions_%s_%s_%f.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), partEta[ieta], suffix.Data()));
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
  histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.51, 49, 1000, 1, 999 );
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{track} (GeV)", "#pi^{#pm}  rejection",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510,42,42);//(#times #epsilon_{pur})
  histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
  histo2DAccEff->GetXaxis()->SetNoExponent();
  histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
  TGraphAsymmErrors* graphPionRejection[maxcalo][nEta]={{NULL}};
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    cout << caloName[icalo].Data() << endl;
    histo2DAccEff->DrawCopy();
    int imatch =1;
    for (Int_t ieta = 0; ieta < nEta; ieta++){
      if(partEta[ieta] > etarangemin[icalo] && partEta[ieta+1] < etarangemax[icalo]){
        cout << partEta[ieta] << " < eta  < " << partEta[ieta+1] << endl;
        if(histPionRejection[imatch][icalo][ieta]){
          graphPionRejection[icalo][ieta] = new TGraphAsymmErrors(histPionRejection[imatch][icalo][ieta]);
          while(graphPionRejection[icalo][ieta]->GetY()[graphPionRejection[icalo][ieta]->GetN()-1] == -1) graphPionRejection[icalo][ieta]->RemovePoint(graphPionRejection[icalo][ieta]->GetN()-1);
          DrawGammaSetMarkerTGraphAsym(graphPionRejection[icalo][ieta], markerStyleEta[ieta], markerSizeEta[ieta], colorEta[ieta] , colorEta[ieta]);
          graphPionRejection[icalo][ieta]->Draw("samelpe");
        }
      }
    }
    TLegend* legendPi0Merge           = GetAndSetLegend2(0.15, 0.92-0.1, 0.3, 0.92-0.1-(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    for (Int_t ieta = 0; ieta < nEta; ieta++){
      if(partEta[ieta] > etarangemin[icalo] && partEta[ieta+1] < etarangemax[icalo]){
        cout << partEta[ieta] << " < eta2  < " << partEta[ieta+1] << endl;
        if(graphPionRejection[icalo][ieta]){
            legendPi0Merge->AddEntry(graphPionRejection[icalo][ieta],Form("%1.1f < #eta < %1.1f", partEta[ieta],partEta[ieta+1]),"p");
        }
      }
    }
    legendPi0Merge->Draw();
    float toplegval = 0.90;
    drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
    drawLatexAdd(collisionSystem.Data(),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
    drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

    // DrawGammaLines(minPtPi0, maxPtPi0woMerged, 1,1,2,kGray+1, 2);


    canvasMergingFraction->Update();
    canvasMergingFraction->Print(Form("%s/EoP_pionrejection_%s.%s",outputDir.Data(),caloName[icalo].Data(),suffix.Data()));
  }
}
