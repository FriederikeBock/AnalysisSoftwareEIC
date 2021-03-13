#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void compare_trackingeffi(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString suffix            = "pdf"
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();


  TString detLabel = "";

//   Double_t maxPtSigma         = 0.175;
  Double_t maxEtaSigma        = 0.005;
  Double_t maxPhiSigma        = 0.01;
  Double_t maxBetaSigma       = 0.05;

  TString outputDir                 = Form("plots/%s/Compare%s",dateForOutput.Data(), addName.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  for (Int_t pid = 0; pid < 6; pid++){
    gSystem->Exec("mkdir -p "+outputDir+"/"+partName[pid]);
  }
  TString collisionSystem = GetCollisionEnergy(addName);
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
  TString outTrackCuts[3]     = {"", "LI2", "LI3"};
  TString labelTrackCuts[3]   = {"", "#geq 2 tracker hits", "#geq 3 tracker hits"};

  TString outBetaCuts[2]      = {"", "LI3"};
  TString labelBetaCuts[2]    = {"", "#geq 3 tracker hits"};
  
  
  TH1D* h_effi_rec_pT[maxNSets][nPID][nEta+1]  = {{{NULL}}};
  TH1D* h_effi_rec_MCpT[maxNSets][nPID][nEta+1] = {{{NULL}}};
  TH1D* h_effi_rec_p[maxNSets][nPID][nEta+1]    = {{{NULL}}};
  TH1D* h_effi_rec_MCp[maxNSets][nPID][nEta+1]   = {{{NULL}}};
  TFile* inputFiles[maxNSets]                 = {NULL};
  TString inputFilesNames[maxNSets];
  TString labels[maxNSets];
  
  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> labels[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< labels[nSets].Data() <<std::endl;
      labels[nSets].ReplaceAll("_"," ");
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;


  for (Int_t iSet = 0; iSet < nSets; iSet++){
    inputFiles[iSet]  = new TFile(inputFilesNames[iSet].Data());
    for (Int_t iEta = 0; iEta < nEta+1; iEta++){
      for (Int_t pid = 0; pid < nPID; pid++){
          h_effi_rec_pT[iSet][pid][iEta]     = (TH1D*)inputFiles[iSet]->Get(Form("TrackingEfficiencyProjections/effi%s_pT_%d", partName[pid].Data(), iEta));
          h_effi_rec_MCpT[iSet][pid][iEta]    = (TH1D*)inputFiles[iSet]->Get(Form("TrackingEfficiencyProjections/effi%s_MCpT_%d", partName[pid].Data(), iEta));
          h_effi_rec_p[iSet][pid][iEta]      = (TH1D*)inputFiles[iSet]->Get(Form("TrackingEfficiencyProjections/effi%s_p_%d", partName[pid].Data(), iEta));
          h_effi_rec_MCp[iSet][pid][iEta]     = (TH1D*)inputFiles[iSet]->Get(Form("TrackingEfficiencyProjections/effi%s_MCp_%d", partName[pid].Data(), iEta));
      }
    }
  }
  
  for(Int_t iEta=0; iEta<nEta+1;iEta++){
    if (!enablePlot[iEta]) continue;
    // 1D PLOT
    Double_t etaMin = partEta[0];
    Double_t etaMax = partEta[nEta];
    if (iEta < nEta){
      etaMin = partEta[iEta];
      etaMax = partEta[iEta+1];
    }

    TCanvas* cReso = new TCanvas("cReso","",0,0,1100,800);
    DrawGammaCanvasSettings( cReso, 0.085, 0.025, 0.02, 0.105);
    // cReso->SetLogz();
  
    TH2F* histoDummyEffiPt   = new TH2F("histoDummyEffiPt","histoDummyEffiPt",1000,0, 20,1000,0.0, 1.25);
    SetStyleHistoTH2ForGraphs(histoDummyEffiPt, "#it{p}_{T} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);
    histoDummyEffiPt->GetXaxis()->SetNoExponent();
    histoDummyEffiPt->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiPt->GetXaxis()->SetMoreLogLabels(kTRUE);
    TH2F* histoDummyEffiMCPt   = new TH2F("histoDummyEffiMCPt","histoDummyEffiMCPt",1000,0, 20,1000,0.0, 1.25);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCPt, "#it{p}_{T}^{MC} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);
    histoDummyEffiMCPt->GetXaxis()->SetNoExponent();
    histoDummyEffiMCPt->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiMCPt->GetXaxis()->SetMoreLogLabels(kTRUE);
    TH2F* histoDummyEffiP   = new TH2F("histoDummyEffiP","histoDummyEffiP",1000,0, 100,1000,0.0, 1.25);
    SetStyleHistoTH2ForGraphs(histoDummyEffiP, "#it{p} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);
    histoDummyEffiP->GetXaxis()->SetNoExponent();
    histoDummyEffiP->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiP->GetXaxis()->SetMoreLogLabels(kTRUE);
    TH2F* histoDummyEffiMCP   = new TH2F("histoDummyEffiMCP","histoDummyEffiMCP",1000,0, 100,1000,0.0, 1.25);
    SetStyleHistoTH2ForGraphs(histoDummyEffiMCP, "#it{p}^{MC} (GeV/#it{c})","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad, 1.15*textSizeSinglePad, 0.9, 0.75);
    histoDummyEffiMCP->GetXaxis()->SetNoExponent();
    histoDummyEffiMCP->GetYaxis()->SetNdivisions(510,kTRUE);
    histoDummyEffiMCP->GetXaxis()->SetMoreLogLabels(kTRUE);
      
    TLegend* legendPtResM  = GetAndSetLegend2(0.4, 0.14, 0.6, 0.14+(nSets*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
    
  
    for (Int_t pid = 0; pid < nPID; pid++){
      histoDummyEffiPt->Draw();
        legendPtResM->Clear();
        for(Int_t iSet=0; iSet<nSets;iSet++){
          DrawGammaSetMarker(h_effi_rec_pT[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_effi_rec_pT[iSet][pid][iEta]->Draw("same,p");
          legendPtResM->AddEntry(h_effi_rec_pT[iSet][pid][iEta],labels[iSet].Data(),"p");
        }
        legendPtResM->Draw();
        drawLatexAdd(collisionSystem,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cReso->Print(Form("%s/%s/Efficiency_pT_%d_%d.%s", outputDir.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  
      
      histoDummyEffiMCPt->Draw();
        for(Int_t iSet=0; iSet<nSets;iSet++){
          DrawGammaSetMarker(h_effi_rec_MCpT[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_effi_rec_MCpT[iSet][pid][iEta]->Draw("same,p");
        }
        legendPtResM->Draw();
        drawLatexAdd(collisionSystem,0.14,0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cReso->Print(Form("%s/%s/Efficiency_MCpT_%d_%d.%s", outputDir.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      
      
      cReso->SetLogx();
      histoDummyEffiP->Draw();
        DrawGammaLines(0, 20, 0., 0., 2, kGray+2, 7);
        for(Int_t iSet=0; iSet<nSets;iSet++){
          DrawGammaSetMarker(h_effi_rec_p[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_effi_rec_p[iSet][pid][iEta]->Draw("same,p");
        }
        legendPtResM->Draw();
        drawLatexAdd(collisionSystem,0.14, 0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cReso->Print(Form("%s/%s/Efficiency_p_%d_%d.%s", outputDir.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
  
      
      histoDummyEffiMCP->Draw();
        for(Int_t iSet=0; iSet<nSets;iSet++){
          DrawGammaSetMarker(h_effi_rec_MCp[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
          h_effi_rec_MCp[iSet][pid][iEta]->Draw("same,p");
        }
        legendPtResM->Draw();
        drawLatexAdd(collisionSystem,0.14, 0.91,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if (pTHard.CompareTo("") != 0) drawLatexAdd(pTHard,0.14, 0.91-textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);    
        drawLatexAdd(Form("%s in %1.1f<#eta<%1.1f", partLabel[pid].Data(), etaMin,etaMax),0.14, 0.91-nLinesCol*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      cReso->Print(Form("%s/%s/Efficiency_MCp_%d_%d.%s", outputDir.Data(), partName[pid].Data(), (Int_t)(etaMin*10), (Int_t)(etaMax*10), suffix.Data()));
      cReso->SetLogx(kFALSE);
    }
  }
}
