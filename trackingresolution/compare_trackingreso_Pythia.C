#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void compare_trackingreso_Pythia(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString suffix            = "pdf",
                            Bool_t properFit          = kTRUE
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();


  TString detLabel = "";

  Double_t maxPtSigma[nEta+1] = {0.42, 0.42, 0.42, 0.42, 0.2, 0.25, 0.2, 0.2, 0.2, 0.25, 
                                 0.15, 0.15, 0.32, 0.42, 0.62, 0.15};

//   Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;
  Double_t maxBetaSigma       = 0.05;
  TString addSigma            = "FitSigma";
  TString addMean             = "FitMean";
  TString addOut              = "Fit";
  if (!properFit){
    maxEtaSigma         = 0.01;
    maxPhiSigma         = 0.1;
    addSigma            = "sigma";
    addMean             = "mean";
    addOut              = "Hist";
    maxBetaSigma        = 0.05;
  }
  TString addLabel        = "";
  // if (primaryTrackSource == 1) {
  //   addLabel += "-innerTracks";
  // }
  TString outputDir                 = Form("plots/%s/Compare%s%s%s",dateForOutput.Data(), addName.Data(), addOut.Data(),addLabel.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  for (Int_t pid = 0; pid < 6; pid++){
    gSystem->Exec("mkdir -p "+outputDir+"/"+partName[pid]);
  }
  TString collisionSystem = "Pythia 6, e+p, 10+250 GeV";
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addName.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }
    
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  Color_t colorSet[8]         = {kGray+1, kBlack, kRed-6, kRed+2, kGreen+1, kCyan+1, kAzure+2, kBlack };
  Style_t markerStyleSet[8]   = {24, 20, 25,  21, 30, 42, 46, 20};
  Size_t markerSizeSet[8]     = {1.5, 1.4, 1.6, 1.5, 1.8, 1.8, 1.5, 1.5 };
  
  const Int_t maxNSets        = 10;

  TString outBetaCuts[2]      = {"", "LI3"};
  TString labelBetaCuts[2]    = {"", "#geq 3 tracker hits"};
  
  TH1D* h_tracks_mean_beta_reso[maxNSets][2][nPID][nEta+1]  = {{{{NULL}}}};
  TH1D* h_tracks_sigma_beta_reso[maxNSets][2][nPID][nEta+1] = {{{{NULL}}}};
  TH1D* h_tracks_mean_pt_reso[maxNSets][nPID][nEta+1]       = {{{NULL}}};
  TH1D* h_tracks_sigma_pt_reso[maxNSets][nPID][nEta+1]      = {{{NULL}}};
  TH1D* h_tracks_mean_p_reso[maxNSets][nPID][nEta+1]         = {{{NULL}}};
  TH1D* h_tracks_sigma_p_reso[maxNSets][nPID][nEta+1]       = {{{NULL}}};
  TH1D* h_tracks_mean_pt_resoEta[maxNSets][nEta+1]          = {{NULL}};
  TH1D* h_tracks_sigma_pt_resoEta[maxNSets][nEta+1]         = {{NULL}};
  TH1D* h_tracks_mean_pt_resoPhi[maxNSets][nEta+1]          = {{NULL}};
  TH1D* h_tracks_sigma_pt_resoPhi[maxNSets][nEta+1]         = {{NULL}};
  TFile* inputFiles[maxNSets]                 = {NULL};
  TString outTrackCutsBase[maxNSets];
  TString inputFilesNames[maxNSets];
  TString labels[maxNSets];
  int primaryTrackSource[maxNSets];
  
  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> labels[nSets] >> primaryTrackSource[nSets] >> outTrackCutsBase[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< labels[nSets].Data()<< "\t"<< primaryTrackSource[nSets] << "\t"<< outTrackCutsBase[nSets].Data() <<  std::endl;
      labels[nSets].ReplaceAll("_"," ");
      if (outTrackCutsBase[nSets] == "All")
        outTrackCutsBase[nSets] = "";
      if(primaryTrackSource[nSets]==1) outTrackCutsBase[nSets]+= "INNER";
      if(primaryTrackSource[nSets]==2) outTrackCutsBase[nSets]+= "SILICON";
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;

  for (Int_t iSet = 0; iSet < nSets; iSet++){
    inputFiles[iSet]  = new TFile(inputFilesNames[iSet].Data());
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      for (Int_t pid = 0; pid < nPID; pid++){
//         cout << Form("histPtResol%s_%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), partName[pid].Data(), addMean.Data() ,iEta,primaryTrackSource[iSet]) << "\t" <<
//         Form("histPtResol%s_%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), partName[pid].Data(), addSigma.Data() ,iEta,primaryTrackSource[iSet]) << endl;
        h_tracks_mean_pt_reso[iSet][pid][iEta]     = (TH1D*)inputFiles[iSet]->Get(Form("histPtResol%s_%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), partName[pid].Data(), addMean.Data() ,iEta,primaryTrackSource[iSet]));
        if(!h_tracks_mean_pt_reso[iSet][pid][iEta]) cout << "not found!" << endl;
        h_tracks_sigma_pt_reso[iSet][pid][iEta]    = (TH1D*)inputFiles[iSet]->Get(Form("histPtResol%s_%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), partName[pid].Data(), addSigma.Data() ,iEta,primaryTrackSource[iSet]));
        h_tracks_mean_p_reso[iSet][pid][iEta]      = (TH1D*)inputFiles[iSet]->Get(Form("histPResol%s_%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), partName[pid].Data(), addMean.Data() ,iEta,primaryTrackSource[iSet]));
        h_tracks_sigma_p_reso[iSet][pid][iEta]     = (TH1D*)inputFiles[iSet]->Get(Form("histPResol%s_%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), partName[pid].Data(), addSigma.Data() ,iEta,primaryTrackSource[iSet]));
        for (Int_t iB = 0; iB< 2; iB++){
          h_tracks_mean_beta_reso[iSet][iB][pid][iEta]   = (TH1D*)inputFiles[iSet]->Get(Form("histBetaResol%s_%s_%s_%d_%d", outBetaCuts[iB].Data(), partName[pid].Data(), addMean.Data() ,iEta,primaryTrackSource[iSet]));
          h_tracks_sigma_beta_reso[iSet][iB][pid][iEta]  = (TH1D*)inputFiles[iSet]->Get(Form("histBetaResol%s_%s_%s_%d_%d", outBetaCuts[iB].Data(), partName[pid].Data(), addSigma.Data() ,iEta,primaryTrackSource[iSet]));
        }
      }
      h_tracks_mean_pt_resoEta[iSet][iEta]    = (TH1D*)inputFiles[iSet]->Get(Form("histEtaResol%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), addMean.Data() ,iEta,primaryTrackSource[iSet]));
      h_tracks_sigma_pt_resoEta[iSet][iEta]   = (TH1D*)inputFiles[iSet]->Get(Form("histEtaResol%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), addSigma.Data() ,iEta,primaryTrackSource[iSet]));
      h_tracks_mean_pt_resoPhi[iSet][iEta]    = (TH1D*)inputFiles[iSet]->Get(Form("histPhiResol%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), addMean.Data() ,iEta,primaryTrackSource[iSet]));
      h_tracks_sigma_pt_resoPhi[iSet][iEta]   = (TH1D*)inputFiles[iSet]->Get(Form("histPhiResol%s_%s_%d_%d", outTrackCutsBase[iSet].Data(), addSigma.Data() ,iEta,primaryTrackSource[iSet]));
    }
  }
  cout << __LINE__ << endl;
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    // 1D PLOT
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }
    cout << "iEta: " << iEta << endl;
    TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
    DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.02, 0.105);
    // cReso->SetLogz();

    TH2F* histoDummyPtResMean   = new TH2F("histoDummyPtResMean","histoDummyPtResMean",1000,0, 20,1000,-0.1, 0.1);
    SetStyleHistoTH2ForGraphs(histoDummyPtResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPtResMean->GetXaxis()->SetNoExponent();
    histoDummyPtResMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPtResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
    TLegend* legendPtResM  = GetAndSetLegend2(0.14, 0.94-(nSets*0.85*textSizeLabelsRel), 0.35, 0.94,0.85*textSizeLabelsPixel, 1, "", 43, 0.1);
    TH2F* histoDummyPtResSigma   = new TH2F("histoDummyPtResSigma","histoDummyPtResSigma",1000,0, 20,1000,-0.0, maxPtSigma[iEta]);
    SetStyleHistoTH2ForGraphs(histoDummyPtResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#it{p}_{T}^{rec} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPtResSigma->GetXaxis()->SetNoExponent();
    histoDummyPtResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPtResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    for (Int_t pid = 0; pid < nPID; pid++){
      histoDummyPtResMean->Draw();
        DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);
        legendPtResM->Clear();
        for(Int_t iSet=0; iSet<nSets;iSet++){
          if (!h_tracks_mean_pt_reso[iSet][pid][iEta]) continue;
          DrawGammaSetMarker(h_tracks_mean_pt_reso[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_tracks_mean_pt_reso[iSet][pid][iEta]->Draw("same,p");
          legendPtResM->AddEntry(h_tracks_mean_pt_reso[iSet][pid][iEta],labels[iSet].Data(),"p");
        }
        legendPtResM->Draw();
        drawLatexAdd(collisionSystem,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.95,0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s/PtResolution_Mean_pT_%d_%d.%s", outputDir.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  
      
      histoDummyPtResSigma->Draw();
        for(Int_t iSet=0; iSet<nSets;iSet++){
          if (!h_tracks_sigma_pt_reso[iSet][pid][iEta]) continue;
          DrawGammaSetMarker(h_tracks_sigma_pt_reso[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_tracks_sigma_pt_reso[iSet][pid][iEta]->Draw("same,p");
        }
        legendPtResM->Draw();
        drawLatexAdd(collisionSystem,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.95,0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s/PtResolution_Sigma_pT_%d_%d.%s", outputDir.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }

    TH2F* histoDummyPResMean   = new TH2F("histoDummyPResMean","histoDummyPResMean",1000,0.1, 220,1000,-0.1, 0.1);
    SetStyleHistoTH2ForGraphs(histoDummyPResMean, "#it{p}^{MC} (GeV/#it{c})","#LT (#it{p}^{rec} - #it{p}^{MC}) / #it{p}^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPResMean->GetXaxis()->SetNoExponent();
    histoDummyPResMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
    TLegend* legendPResM  = GetAndSetLegend2(0.14, 0.94-(nSets*0.85*textSizeLabelsRel), 0.35, 0.94,0.85*textSizeLabelsPixel, 1, "", 43, 0.1);
    TH2F* histoDummyPResSigma   = new TH2F("histoDummyPResSigma","histoDummyPResSigma",1000,0.1, 220,1000,-0.0, maxPtSigma[iEta]);
    SetStyleHistoTH2ForGraphs(histoDummyPResSigma, "#it{p}^{MC} (GeV/#it{c})","#sigma((#it{p}^{rec} - #it{p}^{MC}) / #it{p}^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPResSigma->GetXaxis()->SetNoExponent();
    histoDummyPResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    for (Int_t pid = 0; pid < nPID; pid++){
      cReso->SetLogx();
      histoDummyPResMean->Draw();
        DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);
        legendPResM->Clear();
        for(Int_t iSet=0; iSet<nSets;iSet++){
          if (!h_tracks_mean_p_reso[iSet][pid][iEta]) continue;
          DrawGammaSetMarker(h_tracks_mean_p_reso[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_tracks_mean_p_reso[iSet][pid][iEta]->Draw("same,p");
          legendPResM->AddEntry(h_tracks_mean_p_reso[iSet][pid][iEta],labels[iSet].Data(),"p");
        }
        legendPResM->Draw();
        drawLatexAdd(collisionSystem,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.95,0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s/PResolution_Mean_p_%d_%d.%s", outputDir.Data(), partName[pid].Data(),  (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  
      
      histoDummyPResSigma->Draw();
        for(Int_t iSet=0; iSet<nSets;iSet++){
          if (!h_tracks_sigma_p_reso[iSet][pid][iEta]) continue;
          DrawGammaSetMarker(h_tracks_sigma_p_reso[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_tracks_sigma_p_reso[iSet][pid][iEta]->Draw("same,p");
        }
        legendPResM->Draw();
        drawLatexAdd(collisionSystem,0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.95,0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      cReso->Print(Form("%s/%s/PResolution_Sigma_p_%d_%d.%s", outputDir.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    }
    cReso->SetLogx(kFALSE);
    
    TH2F* histoDummyBetaResMean   = new TH2F("histoDummyBetaResMean","histoDummyBetaResMean",1000,0, 20,1000,-0.01, 0.025);
    SetStyleHistoTH2ForGraphs(histoDummyBetaResMean, "#it{p}^{MC} (GeV/#it{c})","#LT (1/#beta^{rec} - 1/#beta^{MC}) / 1/#beta^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyBetaResMean->GetXaxis()->SetNoExponent();
    histoDummyBetaResMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyBetaResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
    TH2F* histoDummyBetaResSigma   = new TH2F("histoDummyBetaResSigma","histoDummyBetaResSigma",1000,0, 20,1000,-0.0, maxBetaSigma);
    SetStyleHistoTH2ForGraphs(histoDummyBetaResSigma, "#it{p}^{MC} (GeV/#it{c})","#sigma((1/#beta^{rec} - 1/#beta^{MC}) / 1/#beta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyBetaResSigma->GetXaxis()->SetNoExponent();
    histoDummyBetaResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyBetaResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);

    for (Int_t iB = 0; iB < 2; iB++){
      for (Int_t pid = 0; pid < nPID; pid++){
        histoDummyBetaResMean->Draw();
        
          DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);
          legendPtResM->Clear();
          for(Int_t iSet=0; iSet<nSets;iSet++){
            if (h_tracks_mean_beta_reso[iSet][iB][pid][iEta]){
              DrawGammaSetMarker(h_tracks_mean_beta_reso[iSet][iB][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
              h_tracks_mean_beta_reso[iSet][iB][pid][iEta]->Draw("same,p");
              legendPtResM->AddEntry(h_tracks_mean_beta_reso[iSet][iB][pid][iEta],labels[iSet].Data(),"p");
            }
          }
          legendPtResM->Draw();
          drawLatexAdd(collisionSystem,0.95,0.89,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.89-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.95,0.89-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (labelBetaCuts[iB].CompareTo("") != 0) drawLatexAdd( labelBetaCuts[iB].Data(),0.95,0.89-(nLinesCol+1)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        cReso->Print(Form("%s/%s/BetaResolution%s_Mean_pT_%d_%d.%s", outputDir.Data(), partName[pid].Data(), outBetaCuts[iB].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));

        histoDummyBetaResSigma->Draw();
          for(Int_t iSet=0; iSet<nSets;iSet++){
            if (h_tracks_sigma_beta_reso[iSet][iB][pid][iEta]){
              DrawGammaSetMarker(h_tracks_sigma_beta_reso[iSet][iB][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
              h_tracks_sigma_beta_reso[iSet][iB][pid][iEta]->Draw("same,p");
            }
          }
          legendPtResM->Draw();
          drawLatexAdd(collisionSystem,0.95,0.89,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.89-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
          drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.95,0.89-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          if (labelBetaCuts[iB].CompareTo("") != 0) drawLatexAdd( labelBetaCuts[iB].Data(),0.95,0.89-(nLinesCol+1)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);    
        cReso->Print(Form("%s/%s/BetaResolution%s_Sigma_pT_%d_%d.%s", outputDir.Data(), partName[pid].Data(), outBetaCuts[iB].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      }
    }
    
    DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.045, 0.105);
    TH2F* histoDummyEtaResMean   = new TH2F("histoDummyEtaResMean","histoDummyEtaResMean",1000,0, 20,1000,-0.001, 0.001);
    SetStyleHistoTH2ForGraphs(histoDummyEtaResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#eta^{rec} - #eta^{MC}) / #eta^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyEtaResMean->GetXaxis()->SetNoExponent();
    histoDummyEtaResMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEtaResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
    TLegend* legendEtaResM  = GetAndSetLegend2(0.14, 0.91-(nSets*textSizeLabelsRel), 0.55, 0.91,textSizeLabelsPixel, 1, "", 43, 0.1);
    TH2F* histoDummyEtaResSigma   = new TH2F("histoDummyEtaResSigma","histoDummyEtaResSigma",1000,0, 20,1000,-0.0, maxEtaSigma);
    SetStyleHistoTH2ForGraphs(histoDummyEtaResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#eta^{rec} - #eta^{MC}) / #eta^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyEtaResSigma->GetXaxis()->SetNoExponent();
    histoDummyEtaResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEtaResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
      
    
    histoDummyEtaResMean->Draw();
      DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        if (!h_tracks_mean_pt_resoEta[iSet][iEta]) continue;
        DrawGammaSetMarker(h_tracks_mean_pt_resoEta[iSet][iEta],  markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
        h_tracks_mean_pt_resoEta[iSet][iEta]->Draw("same,p");
        legendEtaResM->AddEntry(h_tracks_mean_pt_resoEta[iSet][iEta],labels[iSet].Data(),"p");
      }
      legendEtaResM->Draw();
      drawLatexAdd(collisionSystem,0.95,0.89,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.89-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("(h/e)^{#pm} in %1.1f<#eta<%1.1f", etaMin,etaMax),0.95,0.89-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/EtaResolution_Mean_pT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    histoDummyEtaResSigma->Draw();
      for(Int_t iSet=0; iSet<nSets;iSet++){
        if (!h_tracks_sigma_pt_resoEta[iSet][iEta]) continue;
        DrawGammaSetMarker(h_tracks_sigma_pt_resoEta[iSet][iEta],  markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
        h_tracks_sigma_pt_resoEta[iSet][iEta]->Draw("same,p");
      }
      legendEtaResM->Draw();
      drawLatexAdd(collisionSystem,0.95,0.89,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.89-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("(h/e)^{#pm} in %1.1f<#eta<%1.1f", etaMin,etaMax),0.95,0.89-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/EtaResolution_Sigma_pT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    DrawGammaCanvasSettings( cReso, 0.11, 0.02, 0.045, 0.105);
    TH2F* histoDummyPhiResMean   = new TH2F("histoDummyPhiResMean","histoDummyPhiResMean",1000,0, 20,1000,-0.001, 0.001);
    SetStyleHistoTH2ForGraphs(histoDummyPhiResMean, "#it{p}_{T}^{MC} (GeV/#it{c})","#LT (#varphi^{rec} - #varphi^{MC}) / #varphi^{MC} #GT", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPhiResMean->GetXaxis()->SetNoExponent();
    histoDummyPhiResMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPhiResMean->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoDummyPhiResMean->Draw();
    TLegend* legendPhiResM  = GetAndSetLegend2(0.14, 0.91-(nSets*textSizeLabelsRel), 0.55, 0.91,textSizeLabelsPixel, 1, "", 43, 0.1);
    TH2F* histoDummyPhiResSigma   = new TH2F("histoDummyPhiResSigma","histoDummyPhiResSigma",1000,0, 20,1000,-0.0, maxPhiSigma);
    SetStyleHistoTH2ForGraphs(histoDummyPhiResSigma, "#it{p}_{T}^{MC} (GeV/#it{c})","#sigma((#varphi^{rec} - #varphi^{MC}) / #varphi^{MC})", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.05);
    histoDummyPhiResSigma->GetXaxis()->SetNoExponent();
    histoDummyPhiResSigma->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyPhiResSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoDummyPhiResSigma->Draw();
    
    histoDummyPhiResMean->Draw();
      DrawGammaLines(0, 10, 0., 0., 2, kGray+2, 7);
      for(Int_t iSet=0; iSet<nSets;iSet++){
        if (!h_tracks_mean_pt_resoPhi[iSet][iEta]) continue;
        DrawGammaSetMarker(h_tracks_mean_pt_resoPhi[iSet][iEta],  markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
        h_tracks_mean_pt_resoPhi[iSet][iEta]->Draw("same,p");
        legendPhiResM->AddEntry(h_tracks_mean_pt_resoPhi[iSet][iEta],labels[iSet].Data(),"p");
      }
      legendPhiResM->Draw();
      drawLatexAdd(collisionSystem,0.95,0.89,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.89-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("(h/e)^{#pm} in %1.1f<#eta<%1.1f", etaMin,etaMax),0.95,0.89-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PhiResolution_Mean_pT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
    
    histoDummyPhiResSigma->Draw();
      for(Int_t iSet=0; iSet<nSets;iSet++){
        if (!h_tracks_sigma_pt_resoPhi[iSet][iEta]) continue;
        DrawGammaSetMarker(h_tracks_sigma_pt_resoPhi[iSet][iEta],  markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
        h_tracks_sigma_pt_resoPhi[iSet][iEta]->Draw("same,p");
      }
      legendPhiResM->Draw();
      drawLatexAdd(collisionSystem,0.95,0.89,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.95,0.89-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("(h/e)^{#pm} in %1.1f<#eta<%1.1f", etaMin,etaMax),0.95,0.89-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    cReso->Print(Form("%s/PhiResolution_Sigma_pT_%d_%d.%s", outputDir.Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  }
}
