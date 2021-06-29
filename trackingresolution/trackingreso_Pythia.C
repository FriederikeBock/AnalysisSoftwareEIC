#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))


void trackingreso_Pythia(
                            TString inputFileName   = "file.root",
                            TString suffix          = "pdf",
                            TString addLabel        = "",
                            Bool_t properFit        = kTRUE,
                            Int_t trackCuts         = 0,
                            unsigned short primaryTrackSource = 0
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  SetPlotStyle();
  TString dateForOutput             = ReturnDateStringForOutput();
  TString readTrackClass            = "All";
  TString writeLabel                = "";
  TString labelPlotCuts             = "";
  if (trackCuts == 1){
    readTrackClass                  = "LI2";
    addLabel                        = addLabel+"-LI2";
    writeLabel                      = "LI2";
    labelPlotCuts                   = "#geq 2 tracker hits";
  } else if (trackCuts == 2){
    readTrackClass                  = "LI3";
    addLabel                        = addLabel+"-LI3";
    writeLabel                      = "LI3";
    labelPlotCuts                   = "#geq 3 tracker hits";
  }

  for (Int_t eR = 0; eR < 3; eR++) cout << "before: "<< labelEtaRange[eR].Data() << endl;

  if (primaryTrackSource == 1) {
    addLabel += "-innerTracks";
  }
  
  TString collisionSystem = GetCollisionEnergy(addLabel);
  TString magnetLabel     = GetMagnetLabel(addLabel);
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addLabel.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }
  if (addLabel.Contains("MinQ2_20")) {
    pTHard = "#it{Q}^{2} > 20 GeV^{2}";
    nLinesCol++;
  }
  if (addLabel.Contains("innerTracks")) {
    writeLabel = "INNER";
    labelPlotCuts = "Inner tracks only.";
  }

  TString outputDir                 = Form("plots/%s/%s",dateForOutput.Data(),addLabel.Data());
  TString outputDirPTRes            = Form("plots/%s/%s/PTRes",dateForOutput.Data(),addLabel.Data());
  TString outputDirEtaRes           = Form("plots/%s/%s/EtaRes",dateForOutput.Data(),addLabel.Data());
  TString outputDirPhiRes           = Form("plots/%s/%s/PhiRes",dateForOutput.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDirPhiRes);
  gSystem->Exec("mkdir -p "+outputDirPTRes);
  gSystem->Exec("mkdir -p "+outputDirEtaRes);

  for (Int_t pid = 0; pid < nPID; pid++) gSystem->Exec("mkdir -p "+outputDir+"/"+partName[pid]);
    
  Int_t nActiveEta            = 14;
  Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;
  TString detLabel = GetTrackerLabel(addLabel);
  if (addLabel.Contains("Hist")){
    maxPtSigma         = 0.42;
    maxEtaSigma        = 0.01;
    maxPhiSigma        = 0.1;
  }
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  
  TH1D* h_tracks_mean_pt_reso[nPID][nEta+1]     = {{NULL}};
  TH1D* h_tracks_sigma_pt_reso[nPID][nEta+1]    = {{NULL}};  
  TH1D* h_tracks_mean_p_resoP[nPID][nEta+1]     = {{NULL}};
  TH1D* h_tracks_sigma_p_resoP[nPID][nEta+1]    = {{NULL}};  
  TH1D* h_tracks_mean_pt_resoEta[nEta+1]        = {NULL};
  TH1D* h_tracks_sigma_pt_resoEta[nEta+1]       = {NULL};  
  TH1D* h_tracks_mean_pt_resoPhi[nEta+1]        = {NULL};
  TH1D* h_tracks_sigma_pt_resoPhi[nEta+1]       = {NULL};  

  TH1D* h_tracks_reso_bins[nPID][nEta+1][nPt]   = {{{NULL}}};
  TF1* fit_tracks_reso_bins[nPID][nEta+1][nPt]  = {{{NULL}}};
  TH1D* h_tracks_resoEta_bins[nEta+1][nPt]      = {{NULL}};
  TF1* fit_tracks_resoEta_bins[nEta+1][nPt]     = {{NULL}};
  TH1D* h_tracks_resoPhi_bins[nEta+1][nPt]      = {{NULL}};
  TF1* fit_tracks_resoPhi_bins[nEta+1][nPt]     = {{NULL}};

  TH1D* h_tracks_resoP_bins[nPID][nEta+1][nP]   = {{{NULL}}};
  TF1* fit_tracks_resoP_bins[nPID][nEta+1][nP]  = {{{NULL}}};
  
  TH2F* h_tracks_reso[nPID][nEta+1]             = {{NULL}};
  TH2F* h_tracks_resoP[nPID][nEta+1]            = {{NULL}};
  TH2F* h_tracks_resoEta[nEta+1]                = {NULL};
  TH2F* h_tracks_resoPhi[nEta+1]                = {NULL};
  
  TH2F* h_trackMap_eta_pT[nPID][nCuts]          = {{NULL}};
  TH2F* h_trackMapRatio_eta_pT[nPID][nCuts]     = {{NULL}};
  TH2F* h_trackMapRec_eta_pT[nPID][nCuts]       = {{NULL}};
  TH2F* h_trackMapRecRatio_eta_pT[nPID][nCuts]  = {{NULL}};
  TH2F* h_trackMap_eta_p[nPID][nCuts]           = {{NULL}};
  TH2F* h_trackMapRatio_eta_p[nPID][nCuts]      = {{NULL}};
  TH2F* h_trackMapRec_eta_p[nPID][nCuts]        = {{NULL}};
  TH2F* h_trackMapRecRatio_eta_p[nPID][nCuts]   = {{NULL}};
  
  TFile* inputFile  = new TFile(inputFileName.Data());
  TH1D* histNEvents = (TH1D*)inputFile->Get("nEvents");
  Long_t nEvents    = histNEvents->GetBinContent(1);

  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    for (Int_t pid = 0; pid < nPID; pid++){
      h_tracks_reso[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_tracks_reso_pT_%u_%s_%s_%d", primaryTrackSource, readTrackClass.Data(), partNameET2[pid].Data(), iEta));
      h_tracks_reso[pid][iEta]->Sumw2();
      h_tracks_mean_pt_reso[pid][iEta] = new TH1D(Form("histPtResol_%s_mean_%d", partName[pid].Data(), iEta), 
                                            ";#it{p}_{T}^{MC} (GeV/#it{c}); #LT (#it{p}_{T}^{rec}-#it{p}_{T}^{MC})/#it{p}_{T}^{MC} #GT",
                                            nPt, partPt);
      h_tracks_sigma_pt_reso[pid][iEta] = new TH1D(Form("histPtResol_%s_sigma_%d", partName[pid].Data(), iEta), 
                                              ";#it{p}_{T}^{MC} (GeV/#it{c}); #sigma( (#it{p}_{T}^{rec}-#it{p}_{T}^{MC})/#it{p}_{T}^{MC} )", 
                                              nPt, partPt);

      h_tracks_resoP[pid][iEta]  = (TH2F*)inputFile->Get(Form("h_tracks_reso_p_%u_%s_%s_%d", primaryTrackSource, readTrackClass.Data(), partNameET2[pid].Data(), iEta));
      h_tracks_resoP[pid][iEta]->Sumw2();
      h_tracks_mean_p_resoP[pid][iEta] = new TH1D(Form("histPResol_%s_mean_%d", partName[pid].Data(), iEta), 
                                            ";#it{p}^{MC} (GeV/#it{c}); #LT (#it{p}^{rec}-#it{p}^{MC})/#it{p}^{MC} #GT",
                                            nP, partP);
      h_tracks_sigma_p_resoP[pid][iEta] = new TH1D(Form("histPResol_%s_sigma_%d", partName[pid].Data(), iEta), 
                                              ";#it{p}^{MC} (GeV/#it{c}); #sigma( (#it{p}^{rec}-#it{p}^{MC})/#it{p}^{MC} )", 
                                              nP, partP);

    }
    h_tracks_resoEta[iEta]  = (TH2F*)inputFile->Get(Form("h_tracks_reso_Eta_pT_%u_%s_%d", primaryTrackSource, readTrackClass.Data(), iEta));
    h_tracks_resoEta[iEta]->Sumw2();
    h_tracks_mean_pt_resoEta[iEta] = new TH1D(Form("histPtResolEta_mean_%d", iEta), 
                                           ";#it{p}_{T}^{MC} (GeV/#it{c}); #LT (#eta^{rec}-#eta^{MC})/#eta^{MC} #GT",
                                           nPt, partPt);
    h_tracks_sigma_pt_resoEta[iEta] = new TH1D(Form("histPtResolEta_sigma_%d", iEta), 
                                            ";#it{p}_{T}^{MC} (GeV/#it{c}); #sigma( (#eta^{rec}-#eta^{MC})/#eta^{MC} )", 
                                            nPt, partPt);
    h_tracks_resoPhi[iEta]  = (TH2F*)inputFile->Get(Form("h_tracks_reso_Phi_pT_%u_%s_%d", primaryTrackSource, readTrackClass.Data(), iEta));
    h_tracks_resoPhi[iEta]->Sumw2();
    h_tracks_mean_pt_resoPhi[iEta] = new TH1D(Form("histPtResolPhi_mean_%d", iEta), 
                                           ";#it{p}_{T}^{MC} (GeV/#it{c}); #LT (#varphi^{rec}-#varphi^{MC})/#varphi^{MC} #GT",
                                           nPt, partPt);
    h_tracks_sigma_pt_resoPhi[iEta] = new TH1D(Form("histPtResolPhi_sigma_%d", iEta), 
                                            ";#it{p}_{T}^{MC} (GeV/#it{c}); #sigma( (#varphi^{rec}-#varphi^{MC})/#varphi^{MC} )", 
                                            nPt, partPt);
  }
  
  for (Int_t pid = 0; pid < nPID; pid++){
    for (Int_t cut = 0; cut < nCuts; cut++){
      h_trackMap_eta_pT[pid][cut]   = (TH2F*)inputFile->Get(Form("h_tracks_%u_%s_%s_True_Eta_pT", primaryTrackSource, partNameET2[pid].Data(), nameCuts[cut].Data()));
      h_trackMap_eta_pT[pid][cut]->Sumw2();
      if (cut > 0 ){
        h_trackMapRatio_eta_pT[pid][cut-1] = (TH2F*)h_trackMap_eta_pT[pid][cut]->Clone(Form("h_tracksRatio_%s_%s_True_Eta_pT", partName[pid].Data(), nameCuts[cut].Data()));
        h_trackMapRatio_eta_pT[pid][cut-1]->Sumw2();
        h_trackMapRatio_eta_pT[pid][cut-1]->Divide(h_trackMap_eta_pT[pid][0]);
      }
      if (cut == 11){
        h_trackMapRatio_eta_pT[pid][cut] = (TH2F*)h_trackMap_eta_pT[pid][7]->Clone(Form("h_tracksRatio2_%s_%s_True_Eta_pT", partName[pid].Data(), nameCuts[9].Data()));
        h_trackMapRatio_eta_pT[pid][cut]->Sumw2();
        h_trackMapRatio_eta_pT[pid][cut]->Divide(h_trackMap_eta_pT[pid][9]);
      
      }
      h_trackMapRec_eta_pT[pid][cut]   = (TH2F*)inputFile->Get(Form("h_tracks_%u_%s_%s_Rec_Eta_pT", primaryTrackSource, partNameET2[pid].Data(), nameCuts[cut].Data()));
      h_trackMapRec_eta_pT[pid][cut]->Sumw2();
      if (cut > 0 ){
        h_trackMapRecRatio_eta_pT[pid][cut-1] = (TH2F*)h_trackMapRec_eta_pT[pid][cut]->Clone(Form("h_tracksRatio_%s_%s_Rec_Eta_pT", partName[pid].Data(), nameCuts[cut].Data()));
        h_trackMapRecRatio_eta_pT[pid][cut-1]->Sumw2();
        h_trackMapRecRatio_eta_pT[pid][cut-1]->Divide(h_trackMapRec_eta_pT[pid][0]);
      }
      if (cut == 11){
        h_trackMapRecRatio_eta_pT[pid][cut] = (TH2F*)h_trackMapRec_eta_pT[pid][7]->Clone(Form("h_tracksRatio2_%s_%s_Rec_Eta_pT", partName[pid].Data(), nameCuts[9].Data()));
        h_trackMapRecRatio_eta_pT[pid][cut]->Sumw2();
        h_trackMapRecRatio_eta_pT[pid][cut]->Divide(h_trackMapRec_eta_pT[pid][9]);
      }
      h_trackMap_eta_p[pid][cut]   = (TH2F*)inputFile->Get(Form("h_tracks_%u_%s_%s_True_Eta_p", primaryTrackSource, partNameET2[pid].Data(), nameCuts[cut].Data()));
      h_trackMap_eta_p[pid][cut]->Sumw2();
      if (cut > 0 ){
        h_trackMapRatio_eta_p[pid][cut-1] = (TH2F*)h_trackMap_eta_p[pid][cut]->Clone(Form("h_tracksRatio_%s_%s_True_Eta_p", partName[pid].Data(), nameCuts[cut].Data()));
        h_trackMapRatio_eta_p[pid][cut-1]->Sumw2();
        h_trackMapRatio_eta_p[pid][cut-1]->Divide(h_trackMap_eta_p[pid][0]);
      }
      if (cut == 11){
        h_trackMapRatio_eta_p[pid][cut] = (TH2F*)h_trackMap_eta_p[pid][7]->Clone(Form("h_tracksRatio2_%s_%s_True_Eta_p", partName[pid].Data(), nameCuts[9].Data()));
        h_trackMapRatio_eta_p[pid][cut]->Sumw2();
        h_trackMapRatio_eta_p[pid][cut]->Divide(h_trackMap_eta_p[pid][9]);
        
      }
      h_trackMapRec_eta_p[pid][cut]   = (TH2F*)inputFile->Get(Form("h_tracks_%u_%s_%s_Rec_Eta_p", primaryTrackSource, partNameET2[pid].Data(), nameCuts[cut].Data()));
      h_trackMapRec_eta_p[pid][cut]->Sumw2();
      if (cut > 0 ){
        h_trackMapRecRatio_eta_p[pid][cut-1] = (TH2F*)h_trackMapRec_eta_p[pid][cut]->Clone(Form("h_tracksRatio_%s_%s_Rec_Eta_p", partName[pid].Data(), nameCuts[cut].Data()));
        h_trackMapRecRatio_eta_p[pid][cut-1]->Sumw2();
        h_trackMapRecRatio_eta_p[pid][cut-1]->Divide(h_trackMapRec_eta_p[pid][0]);
      }
      if (cut == 11){
        h_trackMapRecRatio_eta_p[pid][cut] = (TH2F*)h_trackMapRec_eta_p[pid][7]->Clone(Form("h_tracksRatio2_%s_%s_Rec_Eta_p", partName[pid].Data(), nameCuts[9].Data()));
        h_trackMapRecRatio_eta_p[pid][cut]->Sumw2();
        h_trackMapRecRatio_eta_p[pid][cut]->Divide(h_trackMapRec_eta_p[pid][9]);        
      }
    }
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
      for (Int_t pid = 0; pid < nPID; pid++){
  //       std::cout << partPt[iPt] << "\t" << h_tracks_reso[iEta]->GetXaxis()->FindBin(partPt[iPt]) << std::endl;
        h_tracks_reso_bins[pid][iEta][iPt] = (TH1D*)h_tracks_reso[pid][iEta]->ProjectionY(Form("projection_%s_dPt_%d_%d",partName[pid].Data(),iPt, iEta), 
                                                                          h_tracks_reso[pid][iEta]->GetXaxis()->FindBin(partPt[iPt]+0.001), h_tracks_reso[pid][iEta]->GetXaxis()->FindBin(partPt[iPt+1]-0.001),"e"); 
        h_tracks_reso_bins[pid][iEta][iPt]->Rebin(rebinEta_PtResol[iEta]);
        h_tracks_reso_bins[pid][iEta][iPt]->GetXaxis()->SetRangeUser(-0.5,0.5);
        if(h_tracks_reso_bins[pid][iEta][iPt]->GetEntries() > 10 && properFit){
            fit_tracks_reso_bins[pid][iEta][iPt] = new TF1(Form("fitGauss%s%d_%d",partName[pid].Data() ,iPt, iEta),  "gaus", -ptResolFit[iEta],ptResolFit[iEta]);
            fit_tracks_reso_bins[pid][iEta][iPt]->SetParameter(0,h_tracks_reso_bins[pid][iEta][iPt]->GetMaximum() );
            fit_tracks_reso_bins[pid][iEta][iPt]->SetParLimits(0,0.9*h_tracks_reso_bins[pid][iEta][iPt]->GetMaximum(),4*h_tracks_reso_bins[pid][iEta][iPt]->GetMaximum() );
            h_tracks_reso_bins[pid][iEta][iPt]->Fit(fit_tracks_reso_bins[pid][iEta][iPt],"QNRME+","",-ptResolFit[iEta],ptResolFit[iEta]);
        }
      }
      
      h_tracks_resoEta_bins[iEta][iPt] = (TH1D*)h_tracks_resoEta[iEta]->ProjectionY(Form("projection_dEta_%d_%d",iPt, iEta), 
                                                                        h_tracks_resoEta[iEta]->GetXaxis()->FindBin(partPt[iPt]+0.001), h_tracks_resoEta[iEta]->GetXaxis()->FindBin(partPt[iPt+1]-0.001),"e"); 
      h_tracks_resoEta_bins[iEta][iPt]->Rebin(rebinEta_EtaResol[iEta]);
      h_tracks_resoEta_bins[iEta][iPt]->GetXaxis()->SetRangeUser(-0.05,0.05);
      if(h_tracks_resoEta_bins[iEta][iPt]->GetEntries() > 10 && properFit){
          fit_tracks_resoEta_bins[iEta][iPt] = new TF1(Form("fitGaussEta%d_%d",iPt, iEta),  "gaus", -etaResolFit[iEta],etaResolFit[iEta]);
          fit_tracks_resoEta_bins[iEta][iPt]->SetParameter(0,h_tracks_resoEta_bins[iEta][iPt]->GetMaximum() );
          fit_tracks_resoEta_bins[iEta][iPt]->SetParLimits(0,0.9*h_tracks_resoEta_bins[iEta][iPt]->GetMaximum(),4*h_tracks_resoEta_bins[iEta][iPt]->GetMaximum() );
          h_tracks_resoEta_bins[iEta][iPt]->Fit(fit_tracks_resoEta_bins[iEta][iPt],"QNRME+","",-etaResolFit[iEta],etaResolFit[iEta]);
      }
      h_tracks_resoPhi_bins[iEta][iPt] = (TH1D*)h_tracks_resoPhi[iEta]->ProjectionY(Form("projection_dPhi_%d_%d",iPt, iEta), 
                                                                        h_tracks_resoPhi[iEta]->GetXaxis()->FindBin(partPt[iPt]+0.001), h_tracks_resoPhi[iEta]->GetXaxis()->FindBin(partPt[iPt+1]-0.001),"e"); 
      h_tracks_resoPhi_bins[iEta][iPt]->Rebin(rebinEta_PhiResol[iEta]);
      h_tracks_resoPhi_bins[iEta][iPt]->GetXaxis()->SetRangeUser(-0.25,0.25);
      if(h_tracks_resoPhi_bins[iEta][iPt]->GetEntries() > 10 && properFit){
          fit_tracks_resoPhi_bins[iEta][iPt] = new TF1(Form("fitGaussPhi%d_%d",iPt, iEta),  "gaus", -phiResolFit[iEta],phiResolFit[iEta]);
          fit_tracks_resoPhi_bins[iEta][iPt]->SetParameter(0,h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum() );
          fit_tracks_resoPhi_bins[iEta][iPt]->SetParLimits(0,0.9*h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum(),4*h_tracks_resoPhi_bins[iEta][iPt]->GetMaximum() );
          h_tracks_resoPhi_bins[iEta][iPt]->Fit(fit_tracks_resoPhi_bins[iEta][iPt],"QNRME+","",-phiResolFit[iEta],phiResolFit[iEta]);
      }
      
      if(iEta == 0 && debugOutput){
        std::cout <<h_tracks_mean_pt_reso[0][iEta]->GetBinCenter(iPt+1) << "\t"<<iPt+1<<  "\t"<< h_tracks_reso_bins[0][iEta][iPt]->GetMean() <<  "\t"<< h_tracks_reso_bins[0][iEta][iPt]->GetRMS() << std::endl;
        if (fit_tracks_reso_bins[0][iEta][iPt]) std::cout<< fit_tracks_reso_bins[0][iEta][iPt]->GetParameter(1) << "\t" << fit_tracks_reso_bins[0][iEta][iPt]->GetParameter(2) << std::endl;
      }
      if (properFit){
        for (Int_t pid = 0; pid < nPID; pid++){
          if (fit_tracks_reso_bins[pid][iEta][iPt]){
            h_tracks_mean_pt_reso[pid][iEta]->SetBinContent(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParameter(1));
            h_tracks_mean_pt_reso[pid][iEta]->SetBinError(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParError(1));
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinContent(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParameter(2));
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinError(iPt+1,fit_tracks_reso_bins[pid][iEta][iPt]->GetParError(2));
          } else {
            h_tracks_mean_pt_reso[pid][iEta]->SetBinContent(iPt+1,-100000);
            h_tracks_mean_pt_reso[pid][iEta]->SetBinError(iPt+1,0);
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinContent(iPt+1,-1000);
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinError(iPt+1,0);         
          }
        }
        if (fit_tracks_resoEta_bins[iEta][iPt]){
          h_tracks_mean_pt_resoEta[iEta]->SetBinContent(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParameter(1));
          h_tracks_mean_pt_resoEta[iEta]->SetBinError(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParError(1));
          h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParameter(2));
          h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,fit_tracks_resoEta_bins[iEta][iPt]->GetParError(2));          
        } else {
          h_tracks_mean_pt_resoEta[iEta]->SetBinContent(iPt+1,-100000);
          h_tracks_mean_pt_resoEta[iEta]->SetBinError(iPt+1,0);
          h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,-1000);
          h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,0);
        }
        if (fit_tracks_resoPhi_bins[iEta][iPt]){
          h_tracks_mean_pt_resoPhi[iEta]->SetBinContent(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParameter(1));
          h_tracks_mean_pt_resoPhi[iEta]->SetBinError(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParError(1));
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParameter(2));
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,fit_tracks_resoPhi_bins[iEta][iPt]->GetParError(2));          
        } else {
          h_tracks_mean_pt_resoPhi[iEta]->SetBinContent(iPt+1,-100000);
          h_tracks_mean_pt_resoPhi[iEta]->SetBinError(iPt+1,0);
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,-1000);
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,0); 
        }
        
      } else {
        for (Int_t pid = 0; pid < nPID; pid++){
          if(h_tracks_reso_bins[pid][iEta][iPt]->GetEntries() > 10){
            h_tracks_mean_pt_reso[pid][iEta]->SetBinContent(iPt+1,h_tracks_reso_bins[pid][iEta][iPt]->GetMean());
            h_tracks_mean_pt_reso[pid][iEta]->SetBinError(iPt+1,h_tracks_reso_bins[pid][iEta][iPt]->GetMeanError());
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinContent(iPt+1,h_tracks_reso_bins[pid][iEta][iPt]->GetRMS());
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinError(iPt+1,h_tracks_reso_bins[pid][iEta][iPt]->GetRMSError());      
//       h_tracks_sigma_pt_reso[iEta]->SetBinContent(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetStdDev());
//       h_tracks_sigma_pt_reso[iEta]->SetBinError(iPt+1,h_tracks_reso_bins[iEta][iPt]->GetStdDevError());
          } else {
            h_tracks_mean_pt_reso[pid][iEta]->SetBinContent(iPt+1,-100000);
            h_tracks_mean_pt_reso[pid][iEta]->SetBinError(iPt+1,0);
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinContent(iPt+1,-1000);
            h_tracks_sigma_pt_reso[pid][iEta]->SetBinError(iPt+1,0);         
          }
        }        
        if (h_tracks_resoEta_bins[iEta][iPt]->GetEntries() > 10){
          h_tracks_mean_pt_resoEta[iEta]->SetBinContent(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetMean());
          h_tracks_mean_pt_resoEta[iEta]->SetBinError(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetMeanError());
          h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetRMS());
          h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetRMSError());      
//       h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetStdDev());
//       h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,h_tracks_resoEta_bins[iEta][iPt]->GetStdDevError());
        } else {
          h_tracks_mean_pt_resoEta[iEta]->SetBinContent(iPt+1,-100000);
          h_tracks_mean_pt_resoEta[iEta]->SetBinError(iPt+1,0);
          h_tracks_sigma_pt_resoEta[iEta]->SetBinContent(iPt+1,-1000);
          h_tracks_sigma_pt_resoEta[iEta]->SetBinError(iPt+1,0);
        }
        if (h_tracks_resoPhi_bins[iEta][iPt]->GetEntries() > 10){
          h_tracks_mean_pt_resoPhi[iEta]->SetBinContent(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetMean());
          h_tracks_mean_pt_resoPhi[iEta]->SetBinError(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetMeanError());
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetRMS());
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetRMSError());      
//       h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetStdDev());
//       h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,h_tracks_resoPhi_bins[iEta][iPt]->GetStdDevError());
        } else {
          h_tracks_mean_pt_resoPhi[iEta]->SetBinContent(iPt+1,-100000);
          h_tracks_mean_pt_resoPhi[iEta]->SetBinError(iPt+1,0);
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinContent(iPt+1,-1000);
          h_tracks_sigma_pt_resoPhi[iEta]->SetBinError(iPt+1,0);          
        }
      }
    }
    for (Int_t iP = 0; iP < nP; iP++){
      for (Int_t pid = 0; pid < nPID; pid++){
  //       std::cout << partPt[iP] << "\t" << h_tracks_reso[iEta]->GetXaxis()->FindBin(partPt[iP]) << std::endl;
        h_tracks_resoP_bins[pid][iEta][iP] = (TH1D*)h_tracks_resoP[pid][iEta]->ProjectionY(Form("projection_%s_dP_%d_%d",partName[pid].Data(),iP, iEta), 
                                                                          h_tracks_resoP[pid][iEta]->GetXaxis()->FindBin(partP[iP]+0.001), h_tracks_resoP[pid][iEta]->GetXaxis()->FindBin(partP[iP+1]-0.001),"e"); 
        h_tracks_resoP_bins[pid][iEta][iP]->Rebin(rebinEta_PtResol[iEta]);
        h_tracks_resoP_bins[pid][iEta][iP]->GetXaxis()->SetRangeUser(-0.5,0.5);
        if(h_tracks_resoP_bins[pid][iEta][iP]->GetEntries() > 10 && properFit){
            fit_tracks_resoP_bins[pid][iEta][iP] = new TF1(Form("fitGaussP%s%d_%d",partName[pid].Data() ,iP, iEta),  "gaus", -ptResolFit[iEta],ptResolFit[iEta]);
            fit_tracks_resoP_bins[pid][iEta][iP]->SetParameter(0,h_tracks_resoP_bins[pid][iEta][iP]->GetMaximum() );
            fit_tracks_resoP_bins[pid][iEta][iP]->SetParLimits(0,0.9*h_tracks_resoP_bins[pid][iEta][iP]->GetMaximum(),4*h_tracks_resoP_bins[pid][iEta][iP]->GetMaximum() );
            h_tracks_resoP_bins[pid][iEta][iP]->Fit(fit_tracks_resoP_bins[pid][iEta][iP],"QNRME+","",-ptResolFit[iEta],ptResolFit[iEta]);
        }
      }
      if (properFit){
        for (Int_t pid = 0; pid < nPID; pid++){
          if (fit_tracks_resoP_bins[pid][iEta][iP]){
            h_tracks_mean_p_resoP[pid][iEta]->SetBinContent(iP+1,fit_tracks_resoP_bins[pid][iEta][iP]->GetParameter(1));
            h_tracks_mean_p_resoP[pid][iEta]->SetBinError(iP+1,fit_tracks_resoP_bins[pid][iEta][iP]->GetParError(1));
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinContent(iP+1,fit_tracks_resoP_bins[pid][iEta][iP]->GetParameter(2));
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinError(iP+1,fit_tracks_resoP_bins[pid][iEta][iP]->GetParError(2));
          } else {
            h_tracks_mean_p_resoP[pid][iEta]->SetBinContent(iP+1,-100000);
            h_tracks_mean_p_resoP[pid][iEta]->SetBinError(iP+1,0);
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinContent(iP+1,-1000);
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinError(iP+1,0);         
          }
        }
      } else {
        for (Int_t pid = 0; pid < nPID; pid++){
          if(h_tracks_resoP_bins[pid][iEta][iP]->GetEntries() > 10){
            h_tracks_mean_p_resoP[pid][iEta]->SetBinContent(iP+1,h_tracks_resoP_bins[pid][iEta][iP]->GetMean());
            h_tracks_mean_p_resoP[pid][iEta]->SetBinError(iP+1,h_tracks_resoP_bins[pid][iEta][iP]->GetMeanError());
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinContent(iP+1,h_tracks_resoP_bins[pid][iEta][iP]->GetRMS());
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinError(iP+1,h_tracks_resoP_bins[pid][iEta][iP]->GetRMSError());      
//       h_tracks_sigma_pt_reso[iEta]->SetBinContent(iP+1,h_tracks_reso_bins[iEta][iP]->GetStdDev());
//       h_tracks_sigma_pt_reso[iEta]->SetBinError(iP+1,h_tracks_reso_bins[iEta][iP]->GetStdDevError());
          } else {
            h_tracks_mean_p_resoP[pid][iEta]->SetBinContent(iP+1,-100000);
            h_tracks_mean_p_resoP[pid][iEta]->SetBinError(iP+1,0);
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinContent(iP+1,-1000);
            h_tracks_sigma_p_resoP[pid][iEta]->SetBinError(iP+1,0);         
          }
        }        
      }
    }
  }

  if (!trackCuts){
    TCanvas* cSingle2D = new TCanvas("cSingle2D","",0,0,1000,800);
    DrawGammaCanvasSettings( cSingle2D, 0.1, 0.12, 0.025, 0.105);
    cSingle2D->SetLogz();

    for (Int_t pid = 0; pid < nPID; pid++){
      for (Int_t cut = 0; cut < nCuts; cut++){
        cSingle2D->cd();
        cSingle2D->SetLogz();
        SetStyleHistoTH2ForGraphs(h_trackMap_eta_pT[pid][cut], "#it{p}_{T}^{MC} (GeV/#it{c})", "#eta^{MC}", 
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
          h_trackMap_eta_pT[pid][cut]->Scale(1./nEvents);
          h_trackMap_eta_pT[pid][cut]->GetZaxis()->SetRangeUser( 0.9/nEvents, 0.9/nEvents*1e4);
          h_trackMap_eta_pT[pid][cut]->Draw("colz");
          
          drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(2*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(labelCuts[cut].Data(),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
                
        cSingle2D->SaveAs(Form("%s/%s/TrackDistributionPtvsEta_MC_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 

        cSingle2D->cd();
        SetStyleHistoTH2ForGraphs(h_trackMap_eta_p[pid][cut], "#it{p}^{MC} (GeV/#it{c})", "#eta^{MC}", 
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
          h_trackMap_eta_p[pid][cut]->Scale(1./nEvents);
          h_trackMap_eta_p[pid][cut]->GetZaxis()->SetRangeUser( 0.9/nEvents, 0.9/nEvents*1e4);
          h_trackMap_eta_p[pid][cut]->Draw("colz");
          drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(2*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(labelCuts[cut].Data(),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        
        cSingle2D->SaveAs(Form("%s/%s/TrackDistributionPvsEta_MC_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 

        cSingle2D->cd();
        SetStyleHistoTH2ForGraphs(h_trackMapRec_eta_pT[pid][cut], "#it{p}_{T}^{rec} (GeV/#it{c})", "#eta^{rec}", 
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
          h_trackMapRec_eta_pT[pid][cut]->Scale(1./nEvents);
          h_trackMapRec_eta_pT[pid][cut]->Draw("colz");
          h_trackMapRec_eta_pT[pid][cut]->GetZaxis()->SetRangeUser( 0.9/nEvents, 0.9/nEvents*1e4);
          drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(2*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(labelCuts[cut].Data(),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        
        cSingle2D->SaveAs(Form("%s/%s/TrackDistributionPtvsEta_Rec_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 
        
        cSingle2D->cd();
        SetStyleHistoTH2ForGraphs(h_trackMapRec_eta_p[pid][cut], "#it{p}^{rec} (GeV/#it{c})", "#eta^{rec}", 
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
          h_trackMapRec_eta_p[pid][cut]->Scale(1./nEvents);
          h_trackMapRec_eta_p[pid][cut]->GetZaxis()->SetRangeUser( 0.9/nEvents, 0.9/nEvents*1e4);
          h_trackMapRec_eta_p[pid][cut]->Draw("colz");
          
          drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+1)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(2*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(labelCuts[cut].Data(),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

        cSingle2D->SaveAs(Form("%s/%s/TrackDistributionPvsEta_Rec_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 
        
        cSingle2D->cd();
        cSingle2D->SetLogz(0);
        if (cut > 0){
          SetStyleHistoTH2ForGraphs(h_trackMapRatio_eta_pT[pid][cut-1], "#it{p}_{T}^{MC} (GeV/#it{c})", "#eta^{MC}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRatio_eta_pT[pid][cut-1]->Draw("colz");

            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[cut].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[0].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            
    
          cSingle2D->SaveAs(Form("%s/%s/Ratio_TrackDistributionPtvsEta_MC_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 
          cSingle2D->cd();
          SetStyleHistoTH2ForGraphs(h_trackMapRatio_eta_p[pid][cut-1], "#it{p}^{MC} (GeV/#it{c})", "#eta^{MC}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRatio_eta_p[pid][cut-1]->Draw("colz");

            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[cut].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[0].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
          cSingle2D->SaveAs(Form("%s/%s/Ratio_TrackDistributionPvsEta_MC_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 

          
          cSingle2D->cd();
          SetStyleHistoTH2ForGraphs(h_trackMapRecRatio_eta_pT[pid][cut-1], "#it{p}_{T}^{rec} (GeV/#it{c})", "#eta^{rec}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRecRatio_eta_pT[pid][cut-1]->Draw("colz");
            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[cut].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[0].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
          cSingle2D->SaveAs(Form("%s/%s/Ratio_TrackDistributionPtvsEta_Rec_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 

          cSingle2D->cd();
          SetStyleHistoTH2ForGraphs(h_trackMapRecRatio_eta_p[pid][cut-1], "#it{p}^{rec} (GeV/#it{c})", "#eta^{rec}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRecRatio_eta_p[pid][cut-1]->Draw("colz");
            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[cut].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[0].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
          cSingle2D->SaveAs(Form("%s/%s/Ratio_TrackDistributionPvsEta_Rec_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[cut].Data(), suffix.Data())); 
        }
        if (cut == 11){
          SetStyleHistoTH2ForGraphs(h_trackMapRatio_eta_pT[pid][cut], "#it{p}_{T}^{MC} (GeV/#it{c})", "#eta^{MC}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRatio_eta_pT[pid][cut]->Draw("colz");
            
            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[7].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[9].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          cSingle2D->SaveAs(Form("%s/%s/RatioToAnyTime_TrackDistributionPtvsEta_MC_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[7].Data(), suffix.Data())); 
          cSingle2D->cd();
          SetStyleHistoTH2ForGraphs(h_trackMapRatio_eta_p[pid][cut], "#it{p}^{MC} (GeV/#it{c})", "#eta^{MC}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRatio_eta_p[pid][cut]->Draw("colz");
            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[7].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[9].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
          cSingle2D->SaveAs(Form("%s/%s/RatioToAnyTime_TrackDistributionPvsEta_MC_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[7].Data(), suffix.Data())); 

          
          cSingle2D->cd();
          SetStyleHistoTH2ForGraphs(h_trackMapRecRatio_eta_pT[pid][cut], "#it{p}_{T}^{rec} (GeV/#it{c})", "#eta^{rec}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRecRatio_eta_pT[pid][cut]->Draw("colz");
            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[7].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[9].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          cSingle2D->SaveAs(Form("%s/%s/RatioToAnyTime_TrackDistributionPtvsEta_Rec_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[7].Data(), suffix.Data())); 

          cSingle2D->cd();
          SetStyleHistoTH2ForGraphs(h_trackMapRecRatio_eta_p[pid][cut], "#it{p}^{rec} (GeV/#it{c})", "#eta^{rec}", 
                                    0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.93);
            
            h_trackMapRecRatio_eta_p[pid][cut]->Draw("colz");
            drawLatexAdd(collisionSystem,0.85,0.14+(nLinesCol+2)*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.85,0.14+(3*textSizeLabelsRel),textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
            drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.85,0.14+2*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("ratio of %s", labelCuts[7].Data()),0.85,0.14+1*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("to %s", labelCuts[9].Data()),0.85,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
          cSingle2D->SaveAs(Form("%s/%s/RatioToAnyTime_TrackDistributionPvsEta_Rec_%s.%s", outputDir.Data(), partName[pid].Data(), nameCuts[7].Data(), suffix.Data())); 
        }
      }
    }
  }
  
  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  split_canvas(cPNG, "cPNG1", nPt);
  TCanvas* cPNGP; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  split_canvas(cPNGP, "cPNGP", nP);
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    Int_t padnum =0;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    if (debugOutput) std::cout << Form("%1.1f<#eta<%1.1f",etaMin,etaMax)  << std::endl;
      
    for (Int_t pid = 0; pid < nPID; pid++){
      padnum =0;
      for(Int_t iPt=0; iPt< nPt; iPt++){
        cPNG->cd(padnum+1);
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1)->SetLogy(1);
        SetStyleHistoTH1ForGraphs(h_tracks_reso_bins[pid][iEta][iPt], "(#it{p}_{T}^{rec}-#it{p}_{T}^{MC})/#it{p}_{T}^{MC}", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        h_tracks_reso_bins[pid][iEta][iPt]->GetYaxis()->SetRangeUser(0.8,h_tracks_reso_bins[pid][iEta][iPt]->GetMaximum()*2);
        
        DrawGammaSetMarker(h_tracks_reso_bins[pid][iEta][iPt], 20, 1.5, kBlack, kBlack);      
        h_tracks_reso_bins[pid][iEta][iPt]->Draw("p,e");
        
        if (fit_tracks_reso_bins[pid][iEta][iPt] && properFit){
          DrawGammaSetMarkerTF1( fit_tracks_reso_bins[pid][iEta][iPt], 7, 2, kBlue+1);
          fit_tracks_reso_bins[pid][iEta][iPt]->Draw("same");
        }
        
        if (padnum == 0){
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        
        drawLatexAdd(Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f ",partPt[iPt],partPt[iPt+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p}_{T} (GeV/#it{c})<%1.1f \t",partPt[iPt],partPt[iPt+1]) << h_tracks_mean_pt_reso[pid][iEta]->GetBinCenter(iPt+1)<< "\t"<< h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iPt+1) << std::endl;
        DrawGammaLines(h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iPt+1), 
                      h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_tracks_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kGray+2, 2);
        DrawGammaLines(h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_reso[pid][iEta]->GetBinContent(iPt+1), 
                      h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iPt+1)-h_tracks_sigma_pt_reso[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_tracks_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);
        DrawGammaLines(h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_reso[pid][iEta]->GetBinContent(iPt+1), 
                      h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iPt+1)+h_tracks_sigma_pt_reso[pid][iEta]->GetBinContent(iPt+1), 
                      0., 0.05*h_tracks_reso_bins[pid][iEta][iPt]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNG->Print(Form("%s/Tracking_%s_Resolution_%d_%d.%s", outputDirPTRes.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

      cPNGP->cd();
      padnum =0;
      for(Int_t iP=0; iP< nP; iP++){
        cPNGP->cd(padnum+1);
        DrawVirtualPadSettings( cPNGP->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNGP->cd(padnum+1)->SetLogy(1);
        SetStyleHistoTH1ForGraphs(h_tracks_resoP_bins[pid][iEta][iP], "(#it{p}^{rec}-#it{p}^{MC})/#it{p}^{MC}", "counts",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
        h_tracks_resoP_bins[pid][iEta][iP]->GetYaxis()->SetRangeUser(0.8,h_tracks_resoP_bins[pid][iEta][iP]->GetMaximum()*2);
        
        DrawGammaSetMarker(h_tracks_resoP_bins[pid][iEta][iP], 20, 1.5, kBlack, kBlack);      
        h_tracks_resoP_bins[pid][iEta][iP]->Draw("p,e");
        
        if (fit_tracks_resoP_bins[pid][iEta][iP] && properFit){
          DrawGammaSetMarkerTF1( fit_tracks_resoP_bins[pid][iEta][iP], 7, 2, kBlue+1);
          fit_tracks_resoP_bins[pid][iEta][iP]->Draw("same");
        }
        
        if (padnum == 0){
          drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("%1.1f<#eta<%1.1f",etaMin,etaMax),0.16,0.89-1.3*textSizeLabelsRel,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        
        drawLatexAdd(Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f ",partP[iP],partP[iP+1]),0.95,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(iEta == 0 && debugOutput) std::cout << Form("%1.1f<#it{p} (GeV/#it{c})<%1.1f \t",partP[iP],partP[iP+1]) << h_tracks_mean_pt_reso[pid][iEta]->GetBinCenter(iP+1)<< "\t"<< h_tracks_mean_pt_reso[pid][iEta]->GetBinContent(iP+1) << std::endl;
        DrawGammaLines(h_tracks_mean_p_resoP[pid][iEta]->GetBinContent(iP+1), 
                      h_tracks_mean_p_resoP[pid][iEta]->GetBinContent(iP+1), 
                      0., 0.05*h_tracks_resoP_bins[pid][iEta][iP]->GetMaximum(), 1, kGray+2, 2);
        DrawGammaLines(h_tracks_mean_p_resoP[pid][iEta]->GetBinContent(iP+1)-h_tracks_sigma_p_resoP[pid][iEta]->GetBinContent(iP+1), 
                      h_tracks_mean_p_resoP[pid][iEta]->GetBinContent(iP+1)-h_tracks_sigma_p_resoP[pid][iEta]->GetBinContent(iP+1), 
                      0., 0.05*h_tracks_resoP_bins[pid][iEta][iP]->GetMaximum(), 1, kRed+2, 2);
        DrawGammaLines(h_tracks_mean_p_resoP[pid][iEta]->GetBinContent(iP+1)+h_tracks_sigma_p_resoP[pid][iEta]->GetBinContent(iP+1), 
                      h_tracks_mean_p_resoP[pid][iEta]->GetBinContent(iP+1)+h_tracks_sigma_p_resoP[pid][iEta]->GetBinContent(iP+1), 
                      0., 0.05*h_tracks_resoP_bins[pid][iEta][iP]->GetMaximum(), 1, kRed+2, 2);

        padnum++;
      }
      cPNGP->Print(Form("%s/Tracking_%s_ResolutionP_%d_%d.%s", outputDirPTRes.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

    }
    
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
        drawLatexAdd(Form("h^{#pm} in %s",detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
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
    cPNG->Print(Form("%s/Tracking_Resolution_Eta_%d_%d.%s", outputDirEtaRes.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

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
        drawLatexAdd(Form("h^{#pm} in %s",detLabel.Data()),0.16,0.89,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
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
    cPNG->Print(Form("%s/Tracking_Resolution_Phi_%d_%d.%s", outputDirPhiRes.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

  }


  // 1D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
  DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.02, 0.105);
  // cReso->SetLogz();

  TH2F* histoDummyPtResMean   = new TH2F("histoDummyPtResMean","histoDummyPtResMean",1000,0, 20,1000,-0.1, 0.1);
  SetStyleHistoTH2ForGraphs(histoDummyPtResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPtResMean->GetXaxis()->SetNoExponent();
  histoDummyPtResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPtResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyPtResSigma   = new TH2F("histoDummyPtResSigma","histoDummyPtResSigma",1000,0, 20,1000,-0.0, maxPtSigma);
  SetStyleHistoTH2ForGraphs(histoDummyPtResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPtResSigma->GetXaxis()->SetNoExponent();
  histoDummyPtResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPtResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);

  TLegend* legendPtResM  = GetAndSetLegend2(0.14, 0.94-(nActiveEta/3*0.75*textSizeLabelsRel), 0.6, 0.94,0.75*textSizeLabelsPixel, 3, "", 43, 0.2);
  TLegend* legendPtResPID  = GetAndSetLegend2(0.14, 0.94-(3*0.85*textSizeLabelsRel), 0.35, 0.94,0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  
  for (Int_t pid = 0; pid < nPID; pid++){
    histoDummyPtResMean->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_tracks_mean_pt_reso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_tracks_mean_pt_reso[pid][iEta]->Draw("same,p");
      if (pid == 0){
        if (iEta == nEta )
          legendPtResM->AddEntry(h_tracks_mean_pt_reso[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
        else 
          legendPtResM->AddEntry(h_tracks_mean_pt_reso[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
    }
    legendPtResM->Draw();
    DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);
      
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PtResolution_%s_Mean_pT.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));

    histoDummyPtResSigma->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_tracks_sigma_pt_reso[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_tracks_sigma_pt_reso[pid][iEta]->Draw("same,p");
    }
    legendPtResM->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PtResolution_%s_Sigma_pT.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
  }
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    histoDummyPtResSigma->Draw();
    legendPtResPID->Clear();
    for (Int_t pid =0; pid < nPID; pid++){
      
      DrawGammaSetMarker(h_tracks_sigma_pt_reso[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_tracks_sigma_pt_reso[pid][iEta]->Draw("same,p");      
      legendPtResPID->AddEntry(h_tracks_sigma_pt_reso[pid][iEta],partLabel[pid].Data(),"p");
    }
    legendPtResPID->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s", etaMin, etaMax, detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
    cReso->Print(Form("%s/PtResolutionPID_Sigma_pT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  }
  
  TH2F* histoDummyPResMean   = new TH2F("histoDummyPResMean","histoDummyPResMean",1000,0.1, 220,1000,-0.1, 0.1);
  SetStyleHistoTH2ForGraphs(histoDummyPResMean, "#it{p}^{MC} (GeV/#it{c})","#LT (#it{p}^{rec} - #it{p}^{MC}) / #it{p}^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPResMean->GetXaxis()->SetMoreLogLabels();
  histoDummyPResMean->GetXaxis()->SetLabelOffset(0.01);
  histoDummyPResMean->GetXaxis()->SetNoExponent();
  histoDummyPResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  TH2F* histoDummyPResSigma   = new TH2F("histoDummyPResSigma","histoDummyPResSigma",1000,0.1, 220,1000,-0.0, maxPtSigma);
  SetStyleHistoTH2ForGraphs(histoDummyPResSigma, "#it{p}^{MC} (GeV/#it{c})","#sigma((#it{p}^{rec} - #it{p}^{MC}) / #it{p}^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPResSigma->GetXaxis()->SetMoreLogLabels();
  histoDummyPResSigma->GetXaxis()->SetLabelOffset(0.01);
  histoDummyPResSigma->GetXaxis()->SetNoExponent();
  histoDummyPResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);

  for (Int_t eR = 0; eR < 3; eR++) cout << "before: "<< labelEtaRange[eR].Data() << endl;
  
  TLegend* legendPResM  = GetAndSetLegend2(0.14, 0.94-(nActiveEta/3*0.75*textSizeLabelsRel), 0.6, 0.94,0.75*textSizeLabelsPixel, 3, "", 43, 0.2);
  TLegend* legendPResPID  = GetAndSetLegend2(0.14, 0.94-(3*0.85*textSizeLabelsRel), 0.35, 0.94,0.85*textSizeLabelsPixel, 2, "", 43, 0.25);
  TLegend* legendPResEta[3]  = {NULL, NULL, NULL}; 
  cReso->SetLogx();
  for (Int_t pid = 0; pid < nPID; pid++){
    histoDummyPResMean->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_tracks_mean_p_resoP[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_tracks_mean_p_resoP[pid][iEta]->Draw("same,p");
      if (pid == 0){
        if (iEta == nEta )
          legendPResM->AddEntry(h_tracks_mean_p_resoP[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
        else 
          legendPResM->AddEntry(h_tracks_mean_p_resoP[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
      }
    }
    legendPResM->Draw();
    DrawGammaLines(0.08, 220., 0., 0., 2, kGray+2, 7);
      
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
    cReso->Print(Form("%s/PResolution_%s_Mean_p.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));

    histoDummyPResSigma->Draw();
    for(Int_t iEta=0; iEta<nEta+1;iEta++){
      if (!enablePlot[iEta]) continue;
      DrawGammaSetMarker(h_tracks_sigma_p_resoP[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
      h_tracks_sigma_p_resoP[pid][iEta]->Draw("same,p");
    }
    legendPResM->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PResolution_%s_Sigma_p.%s", outputDir.Data(), partName[pid].Data(), suffix.Data()));
    
    for (Int_t eR = 0; eR < 3; eR++){
      histoDummyPResMean->Draw();
      cout << labelEtaRange[eR].Data() << endl;
      if (pid == 0) legendPResEta[eR]  = GetAndSetLegend2(0.14, 0.94-((maxNEtaBins[eR]+1)*0.85*textSizeLabelsRel), 0.35, 0.94,0.85*textSizeLabelsPixel, 1, labelEtaRange[eR], 43, 0.2);
      for(Int_t iEta=minEtaBin[eR]; iEta<maxEtaBin[eR]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_tracks_mean_p_resoP[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_tracks_mean_p_resoP[pid][iEta]->Draw("same,p");
        if (pid == 0){
          if (iEta == nEta )
            legendPResEta[eR]->AddEntry(h_tracks_mean_p_resoP[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[0],partEta[iEta]),"p");
          else 
            legendPResEta[eR]->AddEntry(h_tracks_mean_p_resoP[pid][iEta],Form("%1.1f<#eta<%1.1f",partEta[iEta],partEta[iEta+1]),"p");
        }
      }
      legendPResEta[eR]->Draw();
      DrawGammaLines(0.08, 220., 0., 0., 2, kGray+2, 7);
        
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/PResolution_%s_%s_Mean_p.%s", outputDir.Data(), nameOutEtaRange[eR].Data(), partName[pid].Data(), suffix.Data()));

      histoDummyPResSigma->Draw();
      for(Int_t iEta=minEtaBin[eR]; iEta<maxEtaBin[eR]+1;iEta++){
        if (!enablePlot[iEta]) continue;
        DrawGammaSetMarker(h_tracks_sigma_p_resoP[pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        h_tracks_sigma_p_resoP[pid][iEta]->Draw("same,p");
      }
      legendPResEta[eR]->Draw();
      drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
      drawLatexAdd(Form("%s in %s", partLabel[pid].Data(), detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/PResolution_%s_%s_Sigma_p.%s", outputDir.Data(), nameOutEtaRange[eR].Data() ,partName[pid].Data(), suffix.Data()));
    }
  }
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    histoDummyPResSigma->Draw();
    legendPResPID->Clear();
    for (Int_t pid =0; pid < nPID; pid++){
      
      DrawGammaSetMarker(h_tracks_sigma_p_resoP[pid][iEta], markerStylePID[pid], markerSizePID[pid], colorPID[pid], colorPID[pid]);
      h_tracks_sigma_p_resoP[pid][iEta]->Draw("same,p");      
      legendPResPID->AddEntry(h_tracks_sigma_p_resoP[pid][iEta],partLabel[pid].Data(),"p");
    }
    legendPResPID->Draw();
    drawLatexAdd(collisionSystem,0.95,0.91,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
    drawLatexAdd(Form("%1.1f<#eta<%1.1f, %s", etaMin, etaMax, detLabel.Data()),0.95,0.91-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.91-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.91-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PResolutionPID_Sigma_p_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  }
  cReso->SetLogx(0);
  
  DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.045, 0.105);
  TH2F* histoDummyEtaResMean   = new TH2F("histoDummyEtaResMean","histoDummyEtaResMean",1000,0, 20,1000,-0.001, 0.001);
  SetStyleHistoTH2ForGraphs(histoDummyEtaResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT #eta^{rec} - #eta^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyEtaResMean->GetXaxis()->SetNoExponent();
  histoDummyEtaResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyEtaResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyEtaResMean->Draw();
  TLegend* legendEtaResM  = GetAndSetLegend2(0.14, 0.91-(nActiveEta/3*0.75*textSizeLabelsRel), 0.6, 0.91,0.75*textSizeLabelsPixel, 3, "", 43, 0.2);
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

  drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("(h/e)^{#pm} in %s",detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.88-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/EtaResolution_Mean_pT.%s", outputDir.Data(), suffix.Data()));

  
  TH2F* histoDummyEtaResSigma   = new TH2F("histoDummyEtaResSigma","histoDummyEtaResSigma",1000,0, 20,1000,-0.0, maxEtaSigma);
  SetStyleHistoTH2ForGraphs(histoDummyEtaResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma(#eta^{rec} - #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
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
  drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("(h/e)^{#pm} in %s",detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.88-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/EtaResolution_Sigma_pT.%s", outputDir.Data(), suffix.Data()));

  TH2F* histoDummyPhiResMean   = new TH2F("histoDummyPhiResMean","histoDummyPhiResMean",1000,0, 20,1000,-0.01, 0.01);
  SetStyleHistoTH2ForGraphs(histoDummyPhiResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT #varphi^{rec} - #varphi^{MC} #GT (rad)" , 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
  histoDummyPhiResMean->GetXaxis()->SetNoExponent();
  histoDummyPhiResMean->GetYaxis()->SetNdivisions(510,kTRUE);
  histoDummyPhiResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
  histoDummyPhiResMean->Draw();
  TLegend* legendPhiResM  = GetAndSetLegend2(0.14, 0.91-(nActiveEta/3*0.75*textSizeLabelsRel), 0.6, 0.91,0.75*textSizeLabelsPixel, 3, "", 43, 0.2);
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
  drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("(h/e)^{#pm} in %s",detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.88-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/PhiResolution_Mean_pT.%s", outputDir.Data(), suffix.Data()));

  TH2F* histoDummyPhiResSigma   = new TH2F("histoDummyPhiResSigma","histoDummyPhiResSigma",1000,0, 20,1000,-0.0, maxPhiSigma);
  SetStyleHistoTH2ForGraphs(histoDummyPhiResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma(#varphi^{rec} - #varphi^{MC}) (rad)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
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
  drawLatexAdd(collisionSystem,0.95,0.88,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.88-0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
  drawLatexAdd(Form("(h/e)^{#pm} in %s",detLabel.Data()),0.95,0.88-nLinesCol*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("%s", magnetLabel.Data()),0.95,0.88-(nLinesCol+1)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (writeLabel.CompareTo("") != 0) drawLatexAdd(labelPlotCuts,0.95,0.88-(nLinesCol+2)*0.85*textSizeLabelsRel,0.85*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso->Print(Form("%s/PhiResolution_Sigma_pT.%s", outputDir.Data(), suffix.Data()));
  
  TFile* outputFile  = new TFile(inputFileName.Data(),"UPDATE");
  for (Int_t iEta = 0; iEta < nEta+1; iEta++){
    if (properFit){
      for (Int_t pid = 0; pid < 6; pid++){
        h_tracks_mean_pt_reso[pid][iEta]->Write(Form("histPtResol%s_%s_FitMean_%d", writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
        h_tracks_sigma_pt_reso[pid][iEta]->Write(Form("histPtResol%s_%s_FitSigma_%d",writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
        h_tracks_mean_p_resoP[pid][iEta]->Write(Form("histPResol%s_%s_FitMean_%d", writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
        h_tracks_sigma_p_resoP[pid][iEta]->Write(Form("histPResol%s_%s_FitSigma_%d",writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
      }
      h_tracks_mean_pt_resoEta[iEta]->Write(Form("histEtaResol%s_FitMean_%d", writeLabel.Data(), iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoEta[iEta]->Write(Form("histEtaResol%s_FitSigma_%d",writeLabel.Data(), iEta),TObject::kOverwrite);
      h_tracks_mean_pt_resoPhi[iEta]->Write(Form("histPhiResol%s_FitMean_%d", writeLabel.Data(), iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoPhi[iEta]->Write(Form("histPhiResol%s_FitSigma_%d", writeLabel.Data(), iEta),TObject::kOverwrite);
    } else {
      for (Int_t pid = 0; pid < nPID; pid++){
        h_tracks_mean_pt_reso[pid][iEta]->Write(Form("histPtResol%s_%s_mean_%d", writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
        h_tracks_sigma_pt_reso[pid][iEta]->Write(Form("histPtResol%s_%s_sigma_%d", writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);      
        h_tracks_mean_p_resoP[pid][iEta]->Write(Form("histPResol%s_%s_mean_%d", writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);
        h_tracks_sigma_p_resoP[pid][iEta]->Write(Form("histPResol%s_%s_sigma_%d", writeLabel.Data(), partName[pid].Data(), iEta),TObject::kOverwrite);      
      }
      h_tracks_mean_pt_resoEta[iEta]->Write(Form("histEtaResol%s_mean_%d", writeLabel.Data(), iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoEta[iEta]->Write(Form("histEtaResol%s_sigma_%d", writeLabel.Data(), iEta),TObject::kOverwrite);
      h_tracks_mean_pt_resoPhi[iEta]->Write(Form("histPhiResol%s_mean_%d", writeLabel.Data(), iEta),TObject::kOverwrite);
      h_tracks_sigma_pt_resoPhi[iEta]->Write(Form("histPhiResol%s_sigma_%d", writeLabel.Data(), iEta),TObject::kOverwrite);
    }
  }
  outputFile->Write();
  outputFile->Close();
}
