#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void compare_EnergyResol(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString calo              = "",
                            TString suffix            = "pdf",
                            Int_t exampleBin          = 20
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();
  TString perfLabel                 = "ECCE G4 simulation";

  TString detLabel = "";
  TString outputDir                 = Form("plots/%s/Compare%s",dateForOutput.Data(), addName.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  TString collisionSystem = GetCollisionEnergy(addName);
  TString pTHard = "";
  Int_t nLinesCol = 1;
  if (addName.Contains("pTHard5GeV")) {
    pTHard = "#it{p}_{T}^{hard} #geq 5 GeV/#it{c}";
    nLinesCol++;
  }
  Int_t calDir = 0;
  if (calo.Contains("EEMC") || calo.Contains("EHCAL"))
    calDir = 0;
  else if (calo.Contains("BECAL") || calo.Contains("HCALOUT") || calo.Contains("HCALIN"))
    calDir = 1;
  else if (calo.Contains("FEMC") || calo.Contains("LFHCAL"))  
    calDir = 2;
  
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  Color_t colorSet[8]         = {kGray+1, kBlack, kRed-6, kRed+2, kGreen+1, kCyan+1, kAzure+2, kBlack };
  Style_t markerStyleSet[8]   = {24, 20, 25,  21, 30, 42, 46, 20};
  Size_t markerSizeSet[8]     = {1.5, 1.4, 1.6, 1.5, 1.8, 1.8, 1.5, 1.5 };
  Style_t lineStyleSet[8]     = {1, 7, 8, 9, 4, 5, 3, 6 };
  
  const Int_t maxNSets        = 10;
  const static Double_t partE[]   = {   0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 
                                      5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 
                                      13, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 50, 
                                      67.5, 75, 100, 150
                                  };
  TString readNames[7]              = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_", ""};
  TString labelPart[7]              = {"#gamma", "e^{#pm}", "#pi^{#pm}", "p/#bar{p}", "K^{#pm}", "n", "all"};

  Bool_t have_exampleBin[nPID][nEta+1]                = {{kFALSE}};
  Bool_t have_exampleBinF[nPID]                       = {kFALSE};
  TH1D* h_exampleBin_reso[maxNSets][nPID][nEta+1]     = {{{NULL}}};
  TF1* f_exampleBin_reso[maxNSets][nPID][nEta+1]      = {{{NULL}}};
  TH1D* h_exampleBinF_reso[maxNSets][nPID]            = {{NULL}};
  TF1* f_exampleBinF_reso[maxNSets][nPID]             = {{NULL}};
  TH1D* h_exampleBinFH_reso[maxNSets][nPID]           = {{NULL}};
  TF1* f_exampleBinFH_reso[maxNSets][nPID]            = {{NULL}};
  TFile* inputFiles[maxNSets]                         = {NULL};
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
    for (Int_t pid = 1; pid < nPID; pid++){
      for (Int_t iEta=minEtaBinFull[calDir]; iEta<maxEtaBinFull[calDir]+1;iEta++){
        cout << "trying to find: " <<  Form("projection_dE_%s%d_%d", readNames[pid].Data(), exampleBin, iEta) << endl;
        h_exampleBin_reso[iSet][pid][iEta]      = (TH1D*)inputFiles[iSet]->Get(Form("projection_dE_%s%d_%d", readNames[pid].Data(), exampleBin, iEta));
        if (h_exampleBin_reso[iSet][pid][iEta]){
          have_exampleBin[pid][iEta] = kTRUE;
          cout << "found it" << endl;
        }
        f_exampleBin_reso[iSet][pid][iEta]      = (TF1*)inputFiles[iSet]->Get(Form("fit_reso_%s%d_%d", readNames[pid].Data(), exampleBin, iEta));
      }
      h_exampleBinF_reso[iSet][pid]             = (TH1D*)inputFiles[iSet]->Get(Form("projection_dE_%s%d", readNames[pid].Data(), exampleBin));
      if (h_exampleBinF_reso[iSet][pid]){
        have_exampleBinF[pid] = kTRUE;
      }
      f_exampleBinF_reso[iSet][pid]             = (TF1*)inputFiles[iSet]->Get(Form("fit_reso_%s%d", readNames[pid].Data(), exampleBin));
      h_exampleBinFH_reso[iSet][pid]            = (TH1D*)inputFiles[iSet]->Get(Form("projection_dEhighest_%s%d", readNames[pid].Data(), exampleBin));
      f_exampleBinFH_reso[iSet][pid]            = (TF1*)inputFiles[iSet]->Get(Form("fit_resohighest_%s%d", readNames[pid].Data(), exampleBin));
    }
  }

  
  TCanvas* cExampleBin = new TCanvas("cExampleBin","",0,0,1100,1000);
  DrawGammaCanvasSettings( cExampleBin, 0.12, 0.01, 0.015, 0.105);
  cExampleBin->SetLogy();
  
  for (Int_t pid = 0; pid <nPID-1; pid++){
    for (Int_t iEta=minEtaBinFull[calDir]; iEta<maxEtaBinFull[calDir]+1;iEta++){
      if (!have_exampleBin[pid][iEta]) continue;
      TH1D* histResoEbinshighestFitsExample[maxNSets] = {NULL};
      for (Int_t iSet = 0; iSet < nSets; iSet++){
        if (f_exampleBin_reso[iSet][pid][iEta] )
        histResoEbinshighestFitsExample[iSet]  = (TH1D*)f_exampleBin_reso[iSet][pid][iEta]->GetHistogram();
        for (Int_t k = 0; k < histResoEbinshighestFitsExample[iSet]->GetNbinsX()+1; k++)
          histResoEbinshighestFitsExample[iSet]->SetBinError(k, 0);
      }
    
      for (Int_t iSet = 0; iSet < nSets; iSet++){
        if (histResoEbinshighestFitsExample[iSet])
          histResoEbinshighestFitsExample[iSet]->Scale(1./h_exampleBin_reso[iSet][pid][iEta]->GetEntries());
        h_exampleBin_reso[iSet][pid][iEta]->Scale(1./h_exampleBin_reso[iSet][pid][iEta]->GetEntries());
      }
      
      Float_t maximYRange = 0;
      for (Int_t iSet = 0; iSet < nSets; iSet++){
        if (maximYRange < h_exampleBin_reso[iSet][pid][iEta]->GetMaximum())
          maximYRange = h_exampleBin_reso[iSet][pid][iEta]->GetMaximum();
      }
      
      TLegend* legendExamBin  = GetAndSetLegend2(0.16, 0.88-(textSizeLabelsRel*nSets), 0.5, 0.88,textSizeLabelsRel, 2, "", 42, 0.35);
      
      for (Int_t iSet = 0; iSet < nSets; iSet++){
        SetStyleHistoTH1ForGraphs(h_exampleBin_reso[iSet][pid][iEta], "#it{E}^{rec}/#it{E}^{MC}", "probability",
                                  0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1, 510, 505);
        h_exampleBin_reso[iSet][pid][iEta]->GetXaxis()->SetRangeUser(0.,1.25);
        h_exampleBin_reso[iSet][pid][iEta]->GetYaxis()->SetRangeUser(1e-5*maximYRange,20*maximYRange);
        
        DrawGammaSetMarker(h_exampleBin_reso[iSet][pid][iEta], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
        if (iSet == 0)
          h_exampleBin_reso[iSet][pid][iEta]->Draw("p,e");
        else 
          h_exampleBin_reso[iSet][pid][iEta]->Draw("p,e,same");
        
        if (histResoEbinshighestFitsExample[iSet] ){
          DrawGammaSetMarker(histResoEbinshighestFitsExample[iSet], 0, 0, colorSet[iSet], colorSet[iSet], lineStyleSet[iSet], 3);
          histResoEbinshighestFitsExample[iSet]->Draw("lc,same");
          legendExamBin->AddEntry(histResoEbinshighestFitsExample[iSet], " ", "l");
        }
        legendExamBin->AddEntry(h_exampleBin_reso[iSet][pid][iEta], labels[iSet], "p");
        
      }
      legendExamBin->Draw();
      drawLatexAdd(perfLabel,0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("%s: %1.1f< #eta< %1.1f", calo.Data(), partEta[iEta], partEta[iEta+1]),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#it{E} = %1.1f GeV", partE[exampleBin]),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

      cExampleBin->Print(Form("%s/%s_%s_ExampleBin_%d_Eta%d.%s", outputDir.Data(), calo.Data(), readNames[pid].Data(), exampleBin, iEta, suffix.Data()));    
    }    
    
    if (!have_exampleBinF[pid]) continue;
    TH1D* histResoEbinsFitsExample[maxNSets] = {NULL};
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (f_exampleBinF_reso[iSet][pid] )
      histResoEbinsFitsExample[iSet]  = (TH1D*)f_exampleBinF_reso[iSet][pid]->GetHistogram();
      for (Int_t k = 0; k < histResoEbinsFitsExample[iSet]->GetNbinsX()+1; k++)
        histResoEbinsFitsExample[iSet]->SetBinError(k, 0);
    }
  
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (histResoEbinsFitsExample[iSet])
        histResoEbinsFitsExample[iSet]->Scale(1./h_exampleBinF_reso[iSet][pid]->GetEntries());
      h_exampleBinF_reso[iSet][pid]->Scale(1./h_exampleBinF_reso[iSet][pid]->GetEntries());
    }
    
    Float_t maximYRange = 0;
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (maximYRange < h_exampleBinF_reso[iSet][pid]->GetMaximum())
        maximYRange = h_exampleBinF_reso[iSet][pid]->GetMaximum();
    }
    
    TLegend* legendExamBin2  = GetAndSetLegend2(0.16, 0.88-(textSizeLabelsRel*nSets), 0.5, 0.88,textSizeLabelsRel, 2, "", 42, 0.35);
    
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      SetStyleHistoTH1ForGraphs(h_exampleBinF_reso[iSet][pid], "#it{E}^{rec}/#it{E}^{MC}", "probability",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1, 510, 505);
      h_exampleBinF_reso[iSet][pid]->GetXaxis()->SetRangeUser(0.,1.25);
      h_exampleBinF_reso[iSet][pid]->GetYaxis()->SetRangeUser(1e-5*maximYRange,20*maximYRange);
      
      DrawGammaSetMarker(h_exampleBinF_reso[iSet][pid], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      if (iSet == 0)
        h_exampleBinF_reso[iSet][pid]->Draw("p,e");
      else 
        h_exampleBinF_reso[iSet][pid]->Draw("p,e,same");
      
      if (histResoEbinsFitsExample[iSet] ){
        DrawGammaSetMarker(histResoEbinsFitsExample[iSet], 0, 0, colorSet[iSet], colorSet[iSet], lineStyleSet[iSet], 3);
        histResoEbinsFitsExample[iSet]->Draw("lc,same");
        legendExamBin2->AddEntry(histResoEbinsFitsExample[iSet], " ", "l");
      }
      legendExamBin2->AddEntry(h_exampleBinF_reso[iSet][pid], labels[iSet], "p");
      
    }
    legendExamBin2->Draw();
    drawLatexAdd(perfLabel,0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s: %1.1f< #eta< %1.1f", calo.Data(), nominalEtaRegion[calDir][0], nominalEtaRegion[calDir][1]),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("#it{E} = %1.1f GeV", partE[exampleBin]),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cExampleBin->Print(Form("%s/%s_%s_ExampleBin_%d.%s", outputDir.Data(), calo.Data(), readNames[pid].Data(), exampleBin, suffix.Data()));    

    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (f_exampleBinF_reso[iSet][pid] )
      histResoEbinsFitsExample[iSet]  = (TH1D*)f_exampleBinFH_reso[iSet][pid]->GetHistogram();
      for (Int_t k = 0; k < histResoEbinsFitsExample[iSet]->GetNbinsX()+1; k++)
        histResoEbinsFitsExample[iSet]->SetBinError(k, 0);
    }
  
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (histResoEbinsFitsExample[iSet])
        histResoEbinsFitsExample[iSet]->Scale(1./h_exampleBinFH_reso[iSet][pid]->GetEntries());
      h_exampleBinFH_reso[iSet][pid]->Scale(1./h_exampleBinFH_reso[iSet][pid]->GetEntries());
    }
    
    maximYRange = 0;
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (maximYRange < h_exampleBinFH_reso[iSet][pid]->GetMaximum())
        maximYRange = h_exampleBinFH_reso[iSet][pid]->GetMaximum();
    }
    
    TLegend* legendExamBin3  = GetAndSetLegend2(0.16, 0.88-(textSizeLabelsRel*nSets), 0.5, 0.88,textSizeLabelsRel, 2, "", 42, 0.35);
    
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      SetStyleHistoTH1ForGraphs(h_exampleBinFH_reso[iSet][pid], "#it{E}^{rec}/#it{E}^{MC}", "probability",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1, 510, 505);
      h_exampleBinFH_reso[iSet][pid]->GetXaxis()->SetRangeUser(0.,1.25);
      h_exampleBinFH_reso[iSet][pid]->GetYaxis()->SetRangeUser(1e-5*maximYRange,20*maximYRange);
      
      DrawGammaSetMarker(h_exampleBinFH_reso[iSet][pid], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      if (iSet == 0)
        h_exampleBinFH_reso[iSet][pid]->Draw("p,e");
      else 
        h_exampleBinFH_reso[iSet][pid]->Draw("p,e,same");
      
      if (histResoEbinsFitsExample[iSet] ){
        DrawGammaSetMarker(histResoEbinsFitsExample[iSet], 0, 0, colorSet[iSet], colorSet[iSet], lineStyleSet[iSet], 3);
        histResoEbinsFitsExample[iSet]->Draw("lc,same");
        legendExamBin3->AddEntry(histResoEbinsFitsExample[iSet], " ", "l");
      }
      legendExamBin3->AddEntry(h_exampleBinFH_reso[iSet][pid], labels[iSet], "p");
      
    }
    legendExamBin3->Draw();
    drawLatexAdd(perfLabel,0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s: %1.1f< #eta< %1.1f", calo.Data(), nominalEtaRegion[calDir][0], nominalEtaRegion[calDir][1]),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("#it{E} = %1.1f GeV", partE[exampleBin]),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cExampleBin->Print(Form("%s/%s_%s_ExampleBinHighest_%d.%s", outputDir.Data(), calo.Data(), readNames[pid].Data(), exampleBin, suffix.Data()));    
  }
  
  
  
  
}
