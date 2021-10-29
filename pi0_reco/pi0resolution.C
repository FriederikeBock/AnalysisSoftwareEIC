#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 


double mysine(double* x, double* par) { double 

Amplitude = par[0];
double wavelength = par[1];
double phase = par[2];
return Amplitude*sin(2*TMath::Pi()/wavelength*x[0]+phase); 

}


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
    Bool_t doFitting          = kTRUE
){  

  const int maxcalo = 3;
  TString caloName[maxcalo] = {"FEMC", "EEMC", "BECAL"};
  bool caloactive[maxcalo] = {false};
  Double_t mass_pi0PDG = 0.1349770;

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem           = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString perfLabel                 = "ECCE G4 simulation";
  
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
  const int nEta = 40;
  const int nPhi = 25;
  
  const static Double_t partE[]   = {   0.7, 0.9, 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 4.5,  5,
                                        5.5, 6, 7, 8, 10,
                                       12.5,15,17.5,20,22.5,25,  30,35,  40, 50};
  double energymax = 40.;

  TH2F* histEtrueMinv[maxcalo];
  TH2F* histEtrueMinv_allcls[maxcalo];
  TH2F* histpTMinv;
  TH2F* histErecMinv[maxcalo];
  TH2F* histErecMinv_allcls[maxcalo];

  TH1F* histEtrueMinvbin[maxcalo][50]            = {NULL};
  TH1F* histEtrueMinvbin_allcls[maxcalo][50]            = {NULL};

  TF1*  fithistEtrueMinv[maxcalo][50];
  TF1*  fithistEtrueMinv_allclus[50];
  TF1*  fithistEtrueMinv2[50];
  TF1*  fithistEtrueMinvFINAL[maxcalo][50];

  // TH1F* histResoEVsMinv = new TH1F("histResoEVsMinv","",100,0,41);
  TH1F* histResoEVsMinv[maxcalo];
  // TH1F* histMeanEVsMinv = new TH1F("histMeanEVsMinv","",100,0,41);
  TH1F* histMeanEVsMinv[maxcalo];
  // TH1F* histSigmaEVsMinv = new TH1F("histSigmaEVsMinv","",100,0,41);
  TH1F* histSigmaEVsMinv[maxcalo];
  
  TH1F* histEtrueMinvbinSlice[50]            = {NULL};
  TF1*  fithistEtrueMinvSlice[50];
  TF1*  fithistEtrueMinvSliceFINAL[50];
  
  TH3F* histEtaEMinv;
  TH1F* histEEtaMinvbin[50][50]         = {{NULL}};
  TF1* fitEEtaMinvbin[50][50]           = {{NULL}};
  TF1* fitEEtaMinvbinFINAL[50][50]           = {{NULL}};
  TH1F* histResoEEtaVsMinv[50];
  TH1F* histMeanEEtaVsMinv[50];
  TH1F* histSigmaEEtaVsMinv[50];

  TH2F* histEClusters;
  TH1F* histEClusterbin[50];
  TH1F* histEProb[50];
  TH2F* histM02p1Pt; 
  TH2F* histM02p2Pt; 
  TH2F* histM02SCPt;
  TH2F* histM20p1Pt; 
  TH2F* histM20p2Pt;
  TH2F* histM20SCPt;


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

  int padptbin = 0;
  TF1* fithistPtMinv[50];
  TF1* fithistPtMinvFINAL[50];
  TH1F* histResoPtVsMinv = new TH1F("histResoPtVsMinv","",100,0,31);
  TH1F* histMeanPtVsMinv = new TH1F("histMeanPtVsMinv","",100,0,31);
  TH1F* histSigmaPtVsMinv = new TH1F("histSigmaPtVsMinv","",100,0,31);
  for(int icalo=0;icalo<maxcalo;icalo++){
    cout << "processing " << caloName[icalo].Data() << endl;
    TString NameE_Minv = Form("%s/h_Etrue_Minv_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
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
    TString NameE_Minv_allcls = Form("%s/h_Etrue_Minv_allcls_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    histEtrueMinv_allcls[icalo]    = (TH2F*)inputFile->Get(NameE_Minv_allcls.Data());
    if(!histEtrueMinv_allcls[icalo]){ cout << "\t" <<NameE_Minv_allcls.Data() << " not found" << endl; }
    TString NameErecon_Minv_allcls = Form("%s/h_Erec_Minv_allcls_%s_%s",caloName[icalo].Data(), caloName[icalo].Data(), clusterizerName.Data() );
    histErecMinv_allcls[icalo]    = (TH2F*)inputFile->Get(NameErecon_Minv_allcls.Data());
    if(!histErecMinv_allcls[icalo]){ cout << "\t" <<NameErecon_Minv_allcls.Data() << " not found" << endl; }

    cout << "\tinputs loaded" << endl;

    if (histEtrueMinv[icalo]->GetEntries() > 100){
      TCanvas* Minvpi02D = new TCanvas(Form("Minvpi02D%s",caloName[icalo].Data()),"",0,0,1500,800);
      DrawGammaCanvasSettings( Minvpi02D, 0.095, 0.01, 0.01, 0.105);
      Minvpi02D->Divide(2,1);
      Minvpi02D->cd(1);
      histEtrueMinv[icalo]->GetXaxis()->SetRangeUser(0.,40);
      histEtrueMinv[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      histEtrueMinv_allcls[icalo]->GetXaxis()->SetRangeUser(0.,40);
      histEtrueMinv_allcls[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
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
      Tl1.DrawLatex(0.55*energymax,0.03,Form("Calorimeter: %s ", caloName[icalo].Data()));
      Tl1.DrawLatex(0.55*energymax,0.02,Form("Clusterization type: %s ", clusterizerName.Data()));
      Tl1.DrawLatex(0.55*energymax,0.01,"#pi ^{0} #rightarrow #gamma #gamma");
      Tl1.Draw();
      
      Minvpi02D->Print(Form("%s/Mass_vs_E_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));

      cout << "\tMass_vs_E plot saved" << endl;
      // while(ebin < energymax){
      padebin = 0;

      histResoEVsMinv[icalo] = new TH1F(Form("histResoEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
      histMeanEVsMinv[icalo] = new TH1F(Form("histMeanEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
      histSigmaEVsMinv[icalo] = new TH1F(Form("histSigmaEVsMinv_%s",caloName[icalo].Data()),"",nEne-1,partE);
      for(int ebin=0;ebin<nEne-1;ebin++){
        histEtrueMinvbin[icalo][padebin]    = (TH1F*)histEtrueMinv[icalo]->ProjectionY(Form("projection_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]), 
                                              histEtrueMinv[icalo]->GetXaxis()->FindBin(partE[ebin]), histEtrueMinv[icalo]->GetXaxis()->FindBin(partE[ebin+1]),"e");
        histEtrueMinvbin_allcls[icalo][padebin]    = (TH1F*)histEtrueMinv_allcls[icalo]->ProjectionY(Form("projection_Minv_allcls_%s_%f",caloName[icalo].Data(), partE[ebin]), 
                                              histEtrueMinv_allcls[icalo]->GetXaxis()->FindBin(partE[ebin]), histEtrueMinv_allcls[icalo]->GetXaxis()->FindBin(partE[ebin+1]),"e");

        if(histEtrueMinvbin[icalo][padebin]->GetEntries() > 10){
          histEtrueMinvbin[icalo][padebin]->Rebin(4);
          histEtrueMinvbin_allcls[icalo][padebin]->Rebin(4);
          mean = meanErr = sigma = sigmaErr = 0;

          if (doFitting){   

            fithistEtrueMinv[icalo][padebin]    = new TF1(Form("fit_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]), "crystalball", 0.05, 0.15);
            fithistEtrueMinv[icalo][padebin]->SetParameters(0.55*histEtrueMinvbin[icalo][padebin]->GetMaximum(),histEtrueMinvbin[icalo][padebin]->GetMean(),2.*histEtrueMinvbin[icalo][padebin]->GetRMS(),2.*histEtrueMinvbin[icalo][padebin]->GetRMS(),2.5*histEtrueMinvbin[icalo][padebin]->GetRMS());
            histEtrueMinvbin[icalo][padebin]->Fit(fithistEtrueMinv[icalo][padebin],"L0RMEQ","",minmassfit, maxmassfit);
            
            fithistEtrueMinvFINAL[icalo][padebin] = histEtrueMinvbin[icalo][padebin]->GetFunction(Form("fit_Minv_%s_%f",caloName[icalo].Data(), partE[ebin]));
                
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
      mass2D->Divide(4,5);
      int padnumber = 0;
      TLegend* legend1 = new TLegend(0., 0.01, 0.4, 0.25);
      for(int z=0; z< padebin +1 ; z++ ){
      if(z==3){
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
          //Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2)
          SetStyleHistoTH1ForGraphs( histEtrueMinvbin[icalo][padnumber], "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
          DrawGammaSetMarker(histEtrueMinvbin[icalo][padnumber],34,1, kBlack, kBlack);
          histEtrueMinvbin[icalo][padnumber]->GetYaxis()->SetNdivisions(5);
          histEtrueMinvbin[icalo][padnumber]->GetXaxis()->SetNdivisions(5);
          histEtrueMinvbin[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
          histEtrueMinvbin[icalo][padnumber]->Draw("ep");
          DrawGammaSetMarker(histEtrueMinvbin_allcls[icalo][padnumber],27,1, kGreen+2, kGreen+2);
          histEtrueMinvbin_allcls[icalo][padnumber]->GetXaxis()->SetRangeUser(0,0.3);
          histEtrueMinvbin_allcls[icalo][padnumber]->Draw("same,hist");
          fithistEtrueMinvFINAL[icalo][padnumber]->SetLineColor(kBlue);
          fithistEtrueMinvFINAL[icalo][padnumber]->SetLineWidth(1);
          fithistEtrueMinvFINAL[icalo][padnumber]->Draw("same");
          drawLatexAdd(Form("%1.1f < #it{E} < %1.1f GeV",partE[padnumber],partE[padnumber+1]),0.17,0.87,1.5*textSizeSinglePad,kFALSE);
          DrawGammaLines(allmean[icalo][padnumber], allmean[icalo][padnumber], 0, 0.5*fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(),2, 2);
          DrawGammaLines(allmean[icalo][padnumber] + allsigma[icalo][padnumber], allmean[icalo][padnumber]+ allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
          DrawGammaLines(allmean[icalo][padnumber] - allsigma[icalo][padnumber], allmean[icalo][padnumber]- allsigma[icalo][padnumber], 0, fithistEtrueMinvFINAL[icalo][padnumber]->GetMaximum(), 1, 12);
          if(z==0){
            legend1->AddEntry(histEtrueMinvbin[icalo][padnumber], "Simulation", "p");
            legend1->AddEntry(fithistEtrueMinvFINAL[icalo][padnumber], "Fit", "l");
          }
          padnumber++;
      }
      }

      mass2D->Print(Form("%s/Mass_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));

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
      Tl3.DrawLatex(2,0.17,Form("Calorimeter: %s ", caloName[icalo].Data()));
      Tl3.DrawLatex(2,0.16,Form("Clusterization type: %s ", clusterizerName.Data()));
      Tl3.DrawLatex(2,0.14,"#pi ^{0} #rightarrow #gamma #gamma");
      Tl3.Draw();
      fitresults->cd(2);

      histSigmaEVsMinv[icalo]->GetYaxis()->SetRangeUser(0.,0.2);
      SetStyleHistoTH1ForGraphs(histSigmaEVsMinv[icalo], "E^{MC}", "Width (M_{#gamma #gamma})", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.3);  
      DrawGammaSetMarker(histSigmaEVsMinv[icalo],34,1, kBlack, kBlack);   
      histSigmaEVsMinv[icalo]->GetYaxis()->SetRangeUser(0,0.015); 
      histSigmaEVsMinv[icalo]->Draw("e,p");
      
    
      fitresults->Print(Form("%s/Mean_Width_%s_%s.%s", outputDir.Data(), caloName[icalo].Data(), clusterizerName.Data(), suffix.Data()));
    }
  }

  
  //============= End Fit Results======================================

  float minPtPi0 = 0.6;
  float maxPtPi0woMerged = 39;

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
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin),512,505,42,42);//#it{p}_{T} (GeV/#it{c})
  histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(0.1,19.5);//24.5);
  histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
  histo2DAllPi0FWHM->DrawCopy();

  Color_t colorDet[maxcalo]                          = {kBlack, kRed+2, kBlue+2};//, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2,kOrange+2,kSpring+2};
  Style_t markerStyleDet[maxcalo]   = {20, 47, 34};//, 28, 30, 42, 46, 24, 25, 27, 28, 30};
  Size_t markerSizeDet[maxcalo]     = {2.5, 3, 3};//, 1.5, 1.8, 1.8, 1.5, 1.5, 1.4, 1.9,1.5, 1.8 };
  TGraphAsymmErrors* graphPi0FWHMMeV[maxcalo];
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    if(caloactive[icalo] && histSigmaEVsMinv[icalo]){
      histSigmaEVsMinv[icalo]->Scale(1000);
      graphPi0FWHMMeV[icalo] = new TGraphAsymmErrors(histSigmaEVsMinv[icalo]);
      if(graphPi0FWHMMeV[icalo]){
          DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[icalo], markerStyleDet[icalo], markerSizeDet[icalo], colorDet[icalo] , colorDet[icalo]);
          // DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[0] , colorDet[0]);
          graphPi0FWHMMeV[icalo]->Draw("p,same,e");
          // DrawGammaSetMarkerTGraphAsym(graphPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
          // graphPi0TrueFWHMMeV[i]->Draw("p,same,e");
      }
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
  TLegend* legendMassWidth   = GetAndSetLegend2(0.12, 0.2, 0.32, 0.2+3*textsizeLabelsMass, textSizeLabelsPixel);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    if(caloactive[icalo] && histMeanEVsMinv[icalo]){
      histMeanEVsMinv[icalo]->Scale(1000);
      graphPi0Mass[icalo] = new TGraphAsymmErrors(histMeanEVsMinv[icalo]);
      if(graphPi0Mass[icalo]){
        DrawGammaSetMarkerTGraphAsym(graphPi0Mass[icalo], markerStyleDet[icalo], markerSizeDet[icalo], colorDet[icalo] , colorDet[icalo]);
        graphPi0Mass[icalo]->Draw("p,same,e");
        // DrawGammaSetMarkerTGraphAsym(graphPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
        // graphPi0TrueMass[i]->Draw("p,same,e");
        legendMassWidth->AddEntry(graphPi0Mass[icalo],caloName[icalo].Data(),"p");
      }
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


}
