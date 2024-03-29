#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void trackingreso_flatPt(
                            TString inputFileName   = "file.root",
                            TString suffix          = "pdf",
                            TString addLabel        = "",
                            Bool_t properFit        = kTRUE
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();

  TString outputDir                 = Form("plots/%s/%s",dateForOutput.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TString detLabel = "";
  Bool_t enablePlot[8]        = {1, 1, 1, 1, 1, 1, 1, 1};
  Int_t nActiveEta            = 8;
  Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;
  if (addLabel.Contains("defaultLBL") ){
    detLabel  = "LBL";
//     enablePlot[5] = 0;
    enablePlot[6] = 0;
    nActiveEta    = 8;
  } else if (addLabel.Contains("LBLwithLGAD") ){
    detLabel  = "LBL+FTTL";
    enablePlot[6] = 0;
    nActiveEta    = 8;
  } else if (addLabel.Contains("defaultFST") ){
    detLabel  = "FST";
    enablePlot[6] = 0;
    nActiveEta    = 8;
  } else if (addLabel.Contains("FSTwithLGAD") ){
    detLabel  = "FST+FTTL";
    enablePlot[6] = 0;
    nActiveEta    = 8;
  }
  if (addLabel.Contains("Hist")){
    maxPtSigma         = 0.42;
    maxEtaSigma        = 0.01;
    maxPhiSigma        = 0.1;
  }
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  //************************** Read data **************************************************
  const Int_t nPt                  = 11;
  const static Double_t partPt[]     = {0., 0.5, 1.0,  2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
                                      9.0, 10.0};
//   const Int_t nPt                  = 20;
//   const static Double_t partE[]     = {0., 0.5, 1., 1.5, 2., 2.5, 3.0, 3.5, 4.0, 4.5, 
//                                       5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
//                                       10.};
  const Int_t nEta                  = 7;
  Double_t partEta[8]               = { 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};                                    
  const Int_t rebinEta_PtResol[8]   = { 8, 8, 8, 8, 8, 8, 16,  4};                                    
  const Int_t rebinEta_EtaResol[8]   = { 2, 2, 2, 2, 2, 8, 16,  4};                                    
  const Int_t rebinEta_PhiResol[8]   = { 2, 2, 2, 2, 2, 8, 16,  4};                                    
  Float_t ptResolFit[8]     = { 0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.3, 0.2};                                    
  Float_t etaResolFit[8]    = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.01};                                    
  Float_t phiResolFit[8]    = { 0.05, 0.05, 0.05, 0.05, 0.1, 0.15, 0.2, 0.1};                                    
  
  Color_t colorEta[8]         = {kRed+1, kOrange+7, kOrange, kSpring+5, kGreen+1, kCyan+1, kAzure+2, kBlack };
  Style_t markerStyleEta[8]   = {24, 25, 27, 28, 30, 42, 46, 20};
  Size_t markerSizeEta[8]     = {1.5, 1.4, 1.9, 1.5, 1.8, 1.8, 1.5, 1.5 };
  
  TH1D* h_tracks_mean_pt_reso[nEta+1]       = {NULL};
  TH1D* h_tracks_sigma_pt_reso[nEta+1]      = {NULL};  
  TH1D* h_tracks_mean_pt_resoEta[nEta+1]    = {NULL};
  TH1D* h_tracks_sigma_pt_resoEta[nEta+1]   = {NULL};  
  TH1D* h_tracks_mean_pt_resoPhi[nEta+1]    = {NULL};
  TH1D* h_tracks_sigma_pt_resoPhi[nEta+1]   = {NULL};  
  TH1D* h_tracks_reso_bins[nEta+1][nPt]     = {{NULL}};
  TF1* fit_tracks_reso_bins[nEta+1][nPt]    = {{NULL}};
  TH1D* h_tracks_resoEta_bins[nEta+1][nPt]  = {{NULL}};
  TF1* fit_tracks_resoEta_bins[nEta+1][nPt] = {{NULL}};
  TH1D* h_tracks_resoPhi_bins[nEta+1][nPt]  = {{NULL}};
  TF1* fit_tracks_resoPhi_bins[nEta+1][nPt] = {{NULL}};
  
  TH2F* h_tracks_reso[nEta+1]           = {NULL};
  TH2F* h_tracks_resoEta[nEta+1]        = {NULL};
  TH2F* h_tracks_resoPhi[nEta+1]        = {NULL};
  TH2F* h_tracks_reso_Eta[nPt]          = {NULL};

  TFile* inputFile  = new TFile(inputFileName.Data());
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    h_tracks_reso[iEta]  = (TH2F*)inputFile->Get(Form("h_tracks_reso_pT_%d", iEta));
    h_tracks_reso[iEta]->Sumw2();
    h_tracks_mean_pt_reso[iEta] = new TH1D(Form("histPtResol_mean_%d", iEta), 
                                           ";#it{p}_{T}^{MC} (GeV/#it{c}); #LT (#it{p}_{T}^{rec}-#it{p}_{T}^{MC})/#it{p}_{T}^{MC} #GT",
                                           nPt, partPt);
    h_tracks_sigma_pt_reso[iEta] = new TH1D(Form("histPtResol_sigma_%d", iEta), 
                                            ";#it{p}_{T}^{MC} (GeV/#it{c}); #sigma( (#it{p}_{T}^{rec}-#it{p}_{T}^{MC})/#it{p}_{T}^{MC} )", 
                                            nPt, partPt);
    h_tracks_resoEta[iEta]  = (TH2F*)inputFile->Get(Form("h_tracks_reso_Eta_pT_%d", iEta));
    h_tracks_resoEta[iEta]->Sumw2();
    h_tracks_mean_pt_resoEta[iEta] = new TH1D(Form("histPtResolEta_mean_%d", iEta), 
                                           ";#it{p}_{T}^{MC} (GeV/#it{c}); #LT (#eta^{rec}-#eta^{MC})/#eta^{MC} #GT",
                                           nPt, partPt);
    h_tracks_sigma_pt_resoEta[iEta] = new TH1D(Form("histPtResolEta_sigma_%d", iEta), 
                                            ";#it{p}_{T}^{MC} (GeV/#it{c}); #sigma( (#eta^{rec}-#eta^{MC})/#eta^{MC} )", 
                                            nPt, partPt);
    h_tracks_resoPhi[iEta]  = (TH2F*)inputFile->Get(Form("h_tracks_reso_Phi_pT_%d", iEta));
    h_tracks_resoPhi[iEta]->Sumw2();
    h_tracks_mean_pt_resoPhi[iEta] = new TH1D(Form("histPtResolPhi_mean_%d", iEta), 
                                           ";#it{p}_{T}^{MC} (GeV/#it{c}); #LT (#varphi^{rec}-#varphi^{MC})/#varphi^{MC} #GT",
                                           nPt, partPt);
    h_tracks_sigma_pt_resoPhi[iEta] = new TH1D(Form("histPtResolPhi_sigma_%d", iEta), 
                                            ";#it{p}_{T}^{MC} (GeV/#it{c}); #sigma( (#varphi^{rec}-#varphi^{MC})/#varphi^{MC} )", 
                                            nPt, partPt);
  }
  
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput)std::cout << Form("%1.1f<#eta<%1.1f",etaMin,etaMax)  << std::endl;
    for (Int_t iPt = 0; iPt < nPt; iPt++){
//       std::cout << partPt[iPt] << "\t" << h_tracks_reso[iEta]->GetXaxis()->FindBin(partPt[iPt]) << std::endl;
      h_tracks_reso_bins[iEta][iPt] = (TH1D*)h_tracks_reso[iEta]->ProjectionY(Form("projection_dPt_%d_%d",iPt, iEta), 
                                                                        h_tracks_reso[iEta]->GetXaxis()->FindBin(partPt[iPt]), h_tracks_reso[iEta]->GetXaxis()->FindBin(partPt[iPt+1]),"e"); 
      h_tracks_reso_bins[iEta][iPt]->Rebin(rebinEta_PtResol[iEta]);
      h_tracks_reso_bins[iEta][iPt]->GetXaxis()->SetRangeUser(-0.5,0.5);
      if(h_tracks_reso_bins[iEta][iPt]->GetEntries() > 10 && properFit){
          fit_tracks_reso_bins[iEta][iPt] = new TF1(Form("fitGauss%d_%d",iPt, iEta),  "gaus", -ptResolFit[iEta],ptResolFit[iEta]);
          fit_tracks_reso_bins[iEta][iPt]->SetParameter(0,h_tracks_reso_bins[iEta][iPt]->GetMaximum() );
          fit_tracks_reso_bins[iEta][iPt]->SetParLimits(0,0.9*h_tracks_reso_bins[iEta][iPt]->GetMaximum(),4*h_tracks_reso_bins[iEta][iPt]->GetMaximum() );
          h_tracks_reso_bins[iEta][iPt]->Fit(fit_tracks_reso_bins[iEta][iPt],"QNRME+","",-ptResolFit[iEta],ptResolFit[iEta]);
      }
      h_tracks_resoEta_bins[iEta][iPt] = (TH1D*)h_tracks_resoEta[iEta]->ProjectionY(Form("projection_dEta_%d_%d",iPt, iEta), 
                                                                        h_tracks_resoEta[iEta]->GetXaxis()->FindBin(partPt[iPt]), h_tracks_resoEta[iEta]->GetXaxis()->FindBin(partPt[iPt+1]),"e"); 
      h_tracks_resoEta_bins[iEta][iPt]->Rebin(rebinEta_EtaResol[iEta]);
      h_tracks_resoEta_bins[iEta][iPt]->GetXaxis()->SetRangeUser(-0.05,0.05);
      if(h_tracks_resoEta_bins[iEta][iPt]->GetEntries() > 10 && properFit){
          fit_tracks_resoEta_bins[iEta][iPt] = new TF1(Form("fitGaussEta%d_%d",iPt, iEta),  "gaus", -etaResolFit[iEta],etaResolFit[iEta]);
          fit_tracks_resoEta_bins[iEta][iPt]->SetParameter(0,h_tracks_resoEta_bins[iEta][iPt]->GetMaximum() );
          fit_tracks_resoEta_bins[iEta][iPt]->SetParLimits(0,0.9*h_tracks_resoEta_bins[iEta][iPt]->GetMaximum(),4*h_tracks_resoEta_bins[iEta][iPt]->GetMaximum() );
          h_tracks_resoEta_bins[iEta][iPt]->Fit(fit_tracks_resoEta_bins[iEta][iPt],"QNRME+","",-etaResolFit[iEta],etaResolFit[iEta]);
      }
      h_tracks_resoPhi_bins[iEta][iPt] = (TH1D*)h_tracks_resoPhi[iEta]->ProjectionY(Form("projection_dPhi_%d_%d",iPt, iEta), 
                                                                        h_tracks_resoPhi[iEta]->GetXaxis()->FindBin(partPt[iPt]), h_tracks_resoPhi[iEta]->GetXaxis()->FindBin(partPt[iPt+1]),"e"); 
      h_tracks_resoPhi_bins[iEta][iPt]->Rebin(rebinEta_PhiResol[iEta]);
      h_tracks_resoPhi_bins[iEta][iPt]->GetXaxis()->SetRangeUser(-0.25,0.25);
      if(h_tracks_resoPhi_bins[iEta][iPt]->GetEntries() > 10 && properFit){
          fit_tracks_resoPhi_bins[iEta][iPt] = new TF1(Form("fitGaussPhi%d_%d",iPt, iEta),  "gaus", -phiResolFit[iEta],phiResolFit[iEta]);
          fit_tracks_resoPhi_bins[iEta][iPt]->SetParameter(0,h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum() );
          fit_tracks_resoPhi_bins[iEta][iPt]->SetParLimits(0,0.9*h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum(),4*h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum() );
          h_tracks_resoPhi_bins[iEta][iPt]->Fit(fit_tracks_resoPhi_bins[iEta][iPt],"QNRME+","",-phiResolFit[iEta],phiResolFit[iEta]);
      }
      
      if(iEta == 0 && debugOutput){
        std::cout <<h_tracks_mean_pt_reso[iEta]->GetBinCenter(iPt+1) << "\t"<<iPt+1<<  "\t"<< h_tracks_reso_bins[iEta][iPt]->GetMean() <<  "\t"<< h_tracks_reso_bins[iEta][iPt]->GetRMS() << std::endl;
        if (fit_tracks_reso_bins[iEta][iPt]) std::cout<< fit_tracks_reso_bins[iEta][iPt]->GetParameter(1) << "\t" << fit_tracks_reso_bins[iEta][iPt]->GetParameter(2) << std::endl;
      }
      if (properFit){
        if (fit_tracks_reso_bins[iEta][iPt]){
          h_tracks_mean_pt_reso[iEta]->SetBinContent(iPt+1,fit_tracks_reso_bins[iEta][iPt]->GetParameter(1));
          h_tracks_mean_pt_reso[iEta]->SetBinError(iPt+1,fit_tracks_reso_bins[iEta][iPt]->GetParError(1));
          h_tracks_sigma_pt_reso[iEta]->SetBinContent(iPt+1,fit_tracks_reso_bins[iEta][iPt]->GetParameter(2));
          h_tracks_sigma_pt_reso[iEta]->SetBinError(iPt+1,fit_tracks_reso_bins[iEta][iPt]->GetParError(2));
          h_tracks_mean_pt_resoEta[iEta]->SetBinContent(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParameter(1));
          h_tracks_mean_pt_resoEta[iEta]->SetBinError(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParError(1));
          h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParameter(2));
          h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParError(2));
          h_tracks_mean_pt_resoPhi[iEta]->SetBinContent(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParameter(1));
          h_tracks_mean_pt_resoPhi[iEta]->SetBinError(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParError(1));
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParameter(2));
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParError(2));
        } else {
          h_tracks_mean_pt_reso[iEta]->SetBinContent(iPt+1,-100000);
          h_tracks_mean_pt_reso[iEta]->SetBinError(iPt+1,0);
          h_tracks_sigma_pt_reso[iEta]->SetBinContent(iPt+1,-1000);
          h_tracks_sigma_pt_reso[iEta]->SetBinError(iPt+1,0);
          h_tracks_mean_pt_resoEta[iEta]->SetBinContent(iPt+1,-100000);
          h_tracks_mean_pt_resoEta[iEta]->SetBinError(iPt+1,0);
          h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,-1000);
          h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,0);
          h_tracks_mean_pt_resoPhi[iEta]->SetBinContent(iPt+1,-100000);
          h_tracks_mean_pt_resoPhi[iEta]->SetBinError(iPt+1,0);
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,-1000);
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,0);
        }
      } else {
        h_tracks_mean_pt_reso[iEta]->SetBinContent(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetMean());
        h_tracks_mean_pt_reso[iEta]->SetBinError(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetMeanError());
        h_tracks_sigma_pt_reso[iEta]->SetBinContent(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetRMS());
        h_tracks_sigma_pt_reso[iEta]->SetBinError(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetRMSError());      
//       h_tracks_sigma_pt_reso[iEta]->SetBinContent(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetStdDev());
//       h_tracks_sigma_pt_reso[iEta]->SetBinError(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetStdDevError());
        h_tracks_mean_pt_resoEta[iEta]->SetBinContent(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetMean());
        h_tracks_mean_pt_resoEta[iEta]->SetBinError(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetMeanError());
        h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetRMS());
        h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetRMSError());      
//       h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetStdDev());
//       h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetStdDevError());
        h_tracks_mean_pt_resoPhi[iEta]->SetBinContent(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetMean());
        h_tracks_mean_pt_resoPhi[iEta]->SetBinError(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetMeanError());
        h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetRMS());
        h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetRMSError());      
//       h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetStdDev());
//       h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetStdDevError());
      }
    }
  }

  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  split_canvas(cPNG, "cPNG1", nPt);
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Int_t padnum =0;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput) std::cout << Form("%1.1f<#eta<%1.1f",etaMin,etaMax)  << std::endl;
      
    for(Int_t iPt=0; iPt< nPt; iPt++){
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
      cPNG->cd(padnum+1)->SetLogy(1);
      SetStyleHistoTH1ForGraphs(h_tracks_reso_bins[iEta][iPt], "(#it{p}_{T}^{rec}-#it{p}_{T}^{MC})/#it{p}_{T}^{MC}", "counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
      h_tracks_reso_bins[iEta][iPt]->GetYaxis()->SetRangeUser(0.8,h_tracks_reso_bins[iEta][iPt]->GetMaximum()*2);
      
      DrawGammaSetMarker(h_tracks_reso_bins[iEta][iPt], 20, 1.5, kBlack, kBlack);      
      h_tracks_reso_bins[iEta][iPt]->Draw("p,e");
      
      if (fit_tracks_reso_bins[iEta][iPt] && properFit){
        DrawGammaSetMarkerTF1( fit_tracks_reso_bins[iEta][iPt], 7, 2, kBlue+1);
        fit_tracks_reso_bins[iEta][iPt]->Draw("same");
      }
      
      if (padnum == 0){
        drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
      
      drawLatexAdd(Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f ",partPt[iPt],partPt[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f \t",partPt[iPt],partPt[iPt+1]) << h_tracks_mean_pt_reso[iEta]->GetBinCenter(iPt+1)<< "\t"<< h_tracks_mean_pt_reso[iEta]->GetBinContent(iPt+1) << std::endl;
      DrawGammaLines(h_tracks_mean_pt_reso[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_reso[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_reso_bins[iEta][iPt]->GetMaximum(), 1, kGray+2, 2);
      DrawGammaLines(h_tracks_mean_pt_reso[iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_reso[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_reso[iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_reso[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_reso_bins[iEta][iPt]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(h_tracks_mean_pt_reso[iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_reso[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_reso[iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_reso[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_reso_bins[iEta][iPt]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNG->Print(Form("%s/Tracking_Resolution_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    padnum =0;
    for(Int_t iPt=0; iPt< nPt; iPt++){
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
      cPNG->cd(padnum+1)->SetLogy(1);
      SetStyleHistoTH1ForGraphs(h_tracks_resoEta_bins[iEta][iPt], "(#eta^{rec}-#eta^{MC})/#eta^{MC}","counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
      h_tracks_resoEta_bins[iEta][iPt]->GetYaxis()->SetRangeUser(0.8,h_tracks_resoEta_bins[iEta][iPt]->GetMaximum()*2);
      
      DrawGammaSetMarker(h_tracks_resoEta_bins[iEta][iPt], 20, 1.5, kBlack, kBlack);      
      h_tracks_resoEta_bins[iEta][iPt]->Draw("p,e");
      
      if (fit_tracks_resoEta_bins[iEta][iPt] && properFit){
        DrawGammaSetMarkerTF1( fit_tracks_resoEta_bins[iEta][iPt], 7, 2, kBlue+1);
        fit_tracks_resoEta_bins[iEta][iPt]->Draw("same");
      }
      
      if (padnum == 0){
        drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
      
      drawLatexAdd(Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f ",partPt[iPt],partPt[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f \t",partPt[iPt],partPt[iPt+1]) << h_tracks_mean_pt_resoEta[iEta]->GetBinCenter(iPt+1)<< "\t"<< h_tracks_mean_pt_resoEta[iEta]->GetBinContent(iPt+1) << std::endl;
      DrawGammaLines(h_tracks_mean_pt_resoEta[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_resoEta[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_resoEta_bins[iEta][iPt]->GetMaximum(), 1, kGray+2, 2);
      DrawGammaLines(h_tracks_mean_pt_resoEta[iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_resoEta[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_resoEta[iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_resoEta[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_resoEta_bins[iEta][iPt]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(h_tracks_mean_pt_resoEta[iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_resoEta[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_resoEta[iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_resoEta[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_resoEta_bins[iEta][iPt]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNG->Print(Form("%s/Tracking_Resolution_Eta_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    padnum =0;
    for(Int_t iPt=0; iPt< nPt; iPt++){
      cPNG->cd(padnum+1);
      DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
      cPNG->cd(padnum+1)->SetLogy(1);
      SetStyleHistoTH1ForGraphs(h_tracks_resoPhi_bins[iEta][iPt], "(#varphi^{rec}-#varphi^{MC})/#varphi^{MC}","counts",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
      h_tracks_resoPhi_bins[iEta][iPt]->GetYaxis()->SetRangeUser(0.8,h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum()*2);
      
      DrawGammaSetMarker(h_tracks_resoPhi_bins[iEta][iPt], 20, 1.5, kBlack, kBlack);      
      h_tracks_resoPhi_bins[iEta][iPt]->Draw("p,e");
      
      if (fit_tracks_resoPhi_bins[iEta][iPt] && properFit){
        DrawGammaSetMarkerTF1( fit_tracks_resoPhi_bins[iEta][iPt], 7, 2, kBlue+1);
        fit_tracks_resoPhi_bins[iEta][iPt]->Draw("same");
      }
      
      if (padnum == 0){
        drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
      
      drawLatexAdd(Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f ",partPt[iPt],partPt[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f \t",partPt[iPt],partPt[iPt+1]) << h_tracks_mean_pt_resoPhi[iEta]->GetBinCenter(iPt+1)<< "\t"<< h_tracks_mean_pt_resoPhi[iEta]->GetBinContent(iPt+1) << std::endl;
      DrawGammaLines(h_tracks_mean_pt_resoPhi[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_resoPhi[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum(), 1, kGray+2, 2);
      DrawGammaLines(h_tracks_mean_pt_resoPhi[iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_resoPhi[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_resoPhi[iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_resoPhi[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum(), 1, kRed+2, 2);
      DrawGammaLines(h_tracks_mean_pt_resoPhi[iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_resoPhi[iEta]->GetBinContent(iPt+1), 
                     h_tracks_mean_pt_resoPhi[iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_resoPhi[iEta]->GetBinContent(iPt+1), 
                     0., 0.05*h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum(), 1, kRed+2, 2);

      padnum++;
    }
    cPNG->Print(Form("%s/Tracking_Resolution_Phi_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

  }


  // 1D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
  DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.02, 0.105);
  // cReso->SetLogz();

  TH2F* histoDummyPtResMean   = new TH2F("histoDummyPtResMean","histoDummyPtResMean",1000,0, 10,1000,-0.1, 0.1);
  SetStyleHistoTH2ForGraphs(histoDummyPtResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPtResMean->GetXaxis()->SetNoExponent();
  histoDummyPtResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPtResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyPtResMean->Draw();
  TLegend* legendPtResM  = GetAndSetLegend2(0.14, 0.94-(nActiveEta/2*textSizeLabelsRel), 0.55, 0.94,1.1*textSizeLabelsPixel, 2, "", 43, 0.2);
  DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    DrawGammaSetMarker(h_tracks_mean_pt_reso[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
    h_tracks_mean_pt_reso[iEta]->Draw("same,p");
    if (iEta == nEta)
      legendPtResM->AddEntry(h_tracks_mean_pt_reso[iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
    else 
      legendPtResM->AddEntry(h_tracks_mean_pt_reso[iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
  }
  legendPtResM->Draw();
  drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/PtResolution_Mean_pT.%s", outputDir.Data(), suffix.Data()));

  TH2F* histoDummyPtResSigma   = new TH2F("histoDummyPtResSigma","histoDummyPtResSigma",1000,0, 10,1000,-0.0, maxPtSigma);
  SetStyleHistoTH2ForGraphs(histoDummyPtResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPtResSigma->GetXaxis()->SetNoExponent();
  histoDummyPtResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPtResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyPtResSigma->Draw();
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    DrawGammaSetMarker(h_tracks_sigma_pt_reso[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
    h_tracks_sigma_pt_reso[iEta]->Draw("same,p");
  }
  legendPtResM->Draw();
  drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/PtResolution_Sigma_pT.%s", outputDir.Data(), suffix.Data()));

  DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.045, 0.105);
  TH2F* histoDummyEtaResMean   = new TH2F("histoDummyEtaResMean","histoDummyEtaResMean",1000,0, 10,1000,-0.001, 0.001);
  SetStyleHistoTH2ForGraphs(histoDummyEtaResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#eta^{rec} - #eta^{MC}) / #eta^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyEtaResMean->GetXaxis()->SetNoExponent();
  histoDummyEtaResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEtaResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyEtaResMean->Draw();
  TLegend* legendEtaResM  = GetAndSetLegend2(0.14, 0.91-(nActiveEta/2*textSizeLabelsRel), 0.55, 0.91,1.1*textSizeLabelsPixel, 2, "", 43, 0.2);
  DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    DrawGammaSetMarker(h_tracks_mean_pt_resoEta[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
    h_tracks_mean_pt_resoEta[iEta]->Draw("same,p");
    if (iEta == nEta)
      legendEtaResM->AddEntry(h_tracks_mean_pt_resoEta[iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
    else 
      legendEtaResM->AddEntry(h_tracks_mean_pt_resoEta[iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
  }
  legendEtaResM->Draw();
  drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/EtaResolution_Mean_pT.%s", outputDir.Data(), suffix.Data()));

  
  TH2F* histoDummyEtaResSigma   = new TH2F("histoDummyEtaResSigma","histoDummyEtaResSigma",1000,0, 10,1000,-0.0, maxEtaSigma);
  SetStyleHistoTH2ForGraphs(histoDummyEtaResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#eta^{rec} - #eta^{MC}) / #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyEtaResSigma->GetXaxis()->SetNoExponent();
  histoDummyEtaResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEtaResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyEtaResSigma->Draw();
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    DrawGammaSetMarker(h_tracks_sigma_pt_resoEta[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
    h_tracks_sigma_pt_resoEta[iEta]->Draw("same,p");
  }
  legendEtaResM->Draw();
  drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/EtaResolution_Sigma_pT.%s", outputDir.Data(), suffix.Data()));

  TH2F* histoDummyPhiResMean   = new TH2F("histoDummyPhiResMean","histoDummyPhiResMean",1000,0, 10,1000,-0.01, 0.01);
  SetStyleHistoTH2ForGraphs(histoDummyPhiResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#varphi^{rec} - #varphi^{MC}) / #varphi^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPhiResMean->GetXaxis()->SetNoExponent();
  histoDummyPhiResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPhiResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyPhiResMean->Draw();
  TLegend* legendPhiResM  = GetAndSetLegend2(0.14, 0.91-(nActiveEta/2*textSizeLabelsRel), 0.55, 0.91,1.1*textSizeLabelsPixel, 2, "", 43, 0.2);
  DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    DrawGammaSetMarker(h_tracks_mean_pt_resoPhi[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
    h_tracks_mean_pt_resoPhi[iEta]->Draw("same,p");
    if (iEta == nEta)
      legendPhiResM->AddEntry(h_tracks_mean_pt_resoPhi[iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
    else 
      legendPhiResM->AddEntry(h_tracks_mean_pt_resoPhi[iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
  }
  legendPhiResM->Draw();
  drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/PhiResolution_Mean_pT.%s", outputDir.Data(), suffix.Data()));

  TH2F* histoDummyPhiResSigma   = new TH2F("histoDummyPhiResSigma","histoDummyPhiResSigma",1000,0, 10,1000,-0.0, maxPhiSigma);
  SetStyleHistoTH2ForGraphs(histoDummyPhiResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#varphi^{rec} - #varphi^{MC}) / #varphi^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPhiResSigma->GetXaxis()->SetNoExponent();
  histoDummyPhiResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPhiResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyPhiResSigma->Draw();
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    DrawGammaSetMarker(h_tracks_sigma_pt_resoPhi[iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
    h_tracks_sigma_pt_resoPhi[iEta]->Draw("same,p");
  }
  legendPhiResM->Draw();
  drawLatexAdd(Form("#pi^{-} in %s",detLabel.Data()),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/PhiResolution_Sigma_pT.%s", outputDir.Data(), suffix.Data()));
  
  TFile* outputFile  = new TFile(inputFileName.Data(),"UPDATE");
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    if (properFit){
      h_tracks_mean_pt_reso[iEta]->Write(Form("histPtResol_FitMean_%d", iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_reso[iEta]->Write(Form("histPtResol_FitSigma_%d", iEta),TObject::kOverwrite);
      h_tracks_mean_pt_resoEta[iEta]->Write(Form("histEtaResol_FitMean_%d", iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoEta[iEta]->Write(Form("histEtaResol_FitSigma_%d", iEta),TObject::kOverwrite);
      h_tracks_mean_pt_resoPhi[iEta]->Write(Form("histPhiResol_FitMean_%d", iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoPhi[iEta]->Write(Form("histPhiResol_FitSigma_%d", iEta),TObject::kOverwrite);
    } else {
      h_tracks_mean_pt_reso[iEta]->Write(Form("histPtResol_mean_%d", iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_reso[iEta]->Write(Form("histPtResol_sigma_%d", iEta),TObject::kOverwrite);      
      h_tracks_mean_pt_resoEta[iEta]->Write(Form("histEtaResol_mean_%d", iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoEta[iEta]->Write(Form("histEtaResol_sigma_%d", iEta),TObject::kOverwrite);
      h_tracks_mean_pt_resoPhi[iEta]->Write(Form("histPhiResol_mean_%d", iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoPhi[iEta]->Write(Form("histPhiResol_sigma_%d", iEta),TObject::kOverwrite);
    }
  }
  outputFile->Write();
  outputFile->Close();
}

