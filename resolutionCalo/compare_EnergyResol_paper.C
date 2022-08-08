#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);



void compare_EnergyResol_paper(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString calo              = "",
                            TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();
  TString perfLabel                 = "#it{#bf{ECCE}} simulation";

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
Color_t colorPart[7]              = {807, kBlue+1, kRed+2, kGreen+2, kViolet+2, kGray+2, kBlack};
  Style_t markerStylePart[7]        = {20, 46, 33, 30, 28, 42, 25};
  Style_t lineStylePart[7]          = {3, 4, 5, 7, 9, 3, 1};
  Size_t markerSizePart[7]          = {2.5*1.5, 2*1.5, 3*1.5, 2.7*1.5, 2*1.5, 3*1.5, 2*1.5};

  Bool_t have_exampleBin[nPID][nEta+1]                = {{kFALSE}};
  Bool_t have_exampleBinF[maxNSets][nPID]                       = {{kFALSE}};
  TH1D* h_exampleBin_reso[maxNSets][nPID][nEta+1]     = {{{NULL}}};
  TF1* f_exampleBin_reso[maxNSets][nPID][nEta+1]      = {{{NULL}}};
  TH1D* h_exampleBinF_reso[maxNSets][nPID]            = {{NULL}};
  TF1* f_exampleBinF_reso[maxNSets][nPID]             = {{NULL}};
  TH1D* h_exampleBinFH_reso[maxNSets][nPID]           = {{NULL}};
  TF1* f_exampleBinFH_reso[maxNSets][nPID]            = {{NULL}};
  TFile* inputFiles[maxNSets]                         = {NULL};
  TString inputFilesNames[maxNSets];
  TString labels[maxNSets];
  int exampleBin[maxNSets];
  
  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> labels[nSets] >> exampleBin[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< labels[nSets].Data()  << "\t"<< exampleBin[nSets] <<std::endl;
      labels[nSets].ReplaceAll("_"," ");
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;
  double scaleFac[maxNSets][nPID] = {{1}};
  int nActivePartAll[maxNSets] = {0};
  for (Int_t iSet = 0; iSet < nSets; iSet++){
    inputFiles[iSet]  = new TFile(inputFilesNames[iSet].Data());    
    for (Int_t pid = 0; pid < nPID; pid++){
      // for (Int_t iEta=minEtaBinFull[calDir]; iEta<maxEtaBinFull[calDir]+1;iEta++){
      //   cout << "trying to find: " <<  Form("projection_dE_%s%d_%d", readNames[pid].Data(), exampleBin[iSet], iEta) << endl;
      //   h_exampleBin_reso[iSet][pid][iEta]      = (TH1D*)inputFiles[iSet]->Get(Form("projection_dE_%s%d_%d", readNames[pid].Data(), exampleBin[iSet], iEta));
      //   if (h_exampleBin_reso[iSet][pid][iEta]){
      //     have_exampleBin[pid][iEta] = kTRUE;
      //     cout << "found it" << endl;
      //   }
      //   f_exampleBin_reso[iSet][pid][iEta]      = (TF1*)inputFiles[iSet]->Get(Form("fit_reso_%s%d_%d", readNames[pid].Data(), exampleBin[iSet], iEta));
      // }
      h_exampleBinF_reso[iSet][pid]             = (TH1D*)inputFiles[iSet]->Get(Form("projection_dE_%s%d", readNames[pid].Data(), exampleBin[iSet]));
      if (h_exampleBinF_reso[iSet][pid]){
        have_exampleBinF[iSet][pid] = kTRUE;
        nActivePartAll[iSet]++;
        cout << "found exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << endl;
      }
      f_exampleBinF_reso[iSet][pid]             = (TF1*)inputFiles[iSet]->Get(Form("fit_reso_%s%d", readNames[pid].Data(), exampleBin[iSet]));
      h_exampleBinFH_reso[iSet][pid]            = (TH1D*)inputFiles[iSet]->Get(Form("projection_dEhighest_%s%d", readNames[pid].Data(), exampleBin[iSet]));
      f_exampleBinFH_reso[iSet][pid]            = (TF1*)inputFiles[iSet]->Get(Form("fit_resohighest_%s%d", readNames[pid].Data(), exampleBin[iSet]));
      if(h_exampleBinFH_reso[iSet][pid]){
        scaleFac[iSet][pid] = h_exampleBinFH_reso[iSet][pid]->GetEntries();
      }
    }
  }

  
  TCanvas* cExampleBin = new TCanvas("cExampleBin","",0,0,1100,1000);
  DrawGammaCanvasSettings( cExampleBin, 0.12, 0.01, 0.015, 0.105);
  cExampleBin->SetLogy();
  
    TH1D* histResoEbinsFitsExampleHighest[maxNSets][nPID] = {{NULL}};
  for (Int_t pid = 0; pid <nPID-1; pid++){
    
    TH1D* histResoEbinsFitsExample[maxNSets] = {NULL};
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (!have_exampleBinF[iSet][pid]){
        cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << endl;
        continue;
      }
      if (f_exampleBinF_reso[iSet][pid] ){
        histResoEbinsFitsExample[iSet]  = (TH1D*)f_exampleBinF_reso[iSet][pid]->GetHistogram();
        for (Int_t k = 0; k < histResoEbinsFitsExample[iSet]->GetNbinsX()+1; k++)
          histResoEbinsFitsExample[iSet]->SetBinError(k, 0);

        if (histResoEbinsFitsExample[iSet])
          histResoEbinsFitsExample[iSet]->Scale(1./h_exampleBinF_reso[iSet][pid]->GetEntries());
        h_exampleBinF_reso[iSet][pid]->Scale(1./h_exampleBinF_reso[iSet][pid]->GetEntries());
      }
    }
    
    Float_t maximYRange = 0;
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (!have_exampleBinF[iSet][pid]){
        cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << " for ExampleBin plot" << endl;
        continue;
      }
      if (maximYRange < h_exampleBinF_reso[iSet][pid]->GetMaximum())
        maximYRange = h_exampleBinF_reso[iSet][pid]->GetMaximum();
    }
    
    TLegend* legendExamBin2  = GetAndSetLegend2(0.16, 0.88-(textSizeLabelsRel*nSets), 0.5, 0.88,textSizeLabelsRel, 2, "", 42, 0.35);
    
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (!have_exampleBinF[iSet][pid]){
        cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << endl;
        continue;
      }
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
    // drawLatexAdd(Form("#it{E} = %1.1f GeV", partE[exampleBin[iSet]]),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cExampleBin->Print(Form("%s/%s_%s_ExampleBin.%s", outputDir.Data(), calo.Data(), readNames[pid].Data(), suffix.Data()));    

    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (!have_exampleBinF[iSet][pid]){
        cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << endl;
        continue;
      }
      if (f_exampleBinF_reso[iSet][pid] ){
        histResoEbinsFitsExampleHighest[iSet][pid]  = (TH1D*)f_exampleBinFH_reso[iSet][pid]->GetHistogram();
        for (Int_t k = 0; k < histResoEbinsFitsExampleHighest[iSet][pid]->GetNbinsX()+1; k++)
          histResoEbinsFitsExampleHighest[iSet][pid]->SetBinError(k, 0);
    
        if (histResoEbinsFitsExampleHighest[iSet][pid])
          histResoEbinsFitsExampleHighest[iSet][pid]->Scale(1./scaleFac[iSet][pid]);
        h_exampleBinFH_reso[iSet][pid]->Scale(1./scaleFac[iSet][pid]);
      }
    }
    maximYRange = 0;
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (!have_exampleBinF[iSet][pid]){
        cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << endl;
        continue;
      }
      if (maximYRange < h_exampleBinFH_reso[iSet][pid]->GetMaximum())
        maximYRange = h_exampleBinFH_reso[iSet][pid]->GetMaximum();
    }
    
    TLegend* legendExamBin3  = GetAndSetLegend2(0.16, 0.88-(textSizeLabelsRel*nSets), 0.5, 0.88,textSizeLabelsRel, 2, "", 42, 0.35);
    
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if (!have_exampleBinF[iSet][pid]){
        cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << " for ExampleBinHighest plot"  << endl;
        continue;
      }
      SetStyleHistoTH1ForGraphs(h_exampleBinFH_reso[iSet][pid], "#it{E}^{rec}/#it{E}^{MC}", "probability",
                                0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1, 510, 505);
      h_exampleBinFH_reso[iSet][pid]->GetXaxis()->SetRangeUser(0.,1.25);
      h_exampleBinFH_reso[iSet][pid]->GetYaxis()->SetRangeUser(1e-5*maximYRange,20*maximYRange);
      
      DrawGammaSetMarker(h_exampleBinFH_reso[iSet][pid], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      if (iSet == 0)
        h_exampleBinFH_reso[iSet][pid]->Draw("p,e");
      else 
        h_exampleBinFH_reso[iSet][pid]->Draw("p,e,same");
      
      if (histResoEbinsFitsExampleHighest[iSet][pid] ){
        DrawGammaSetMarker(histResoEbinsFitsExampleHighest[iSet][pid], 0, 0, colorSet[iSet], colorSet[iSet], lineStyleSet[iSet], 3);
        histResoEbinsFitsExampleHighest[iSet][pid]->Draw("lc,same");
        legendExamBin3->AddEntry(histResoEbinsFitsExampleHighest[iSet][pid], " ", "l");
      }
      legendExamBin3->AddEntry(h_exampleBinFH_reso[iSet][pid], labels[iSet], "p");
      
    }
    legendExamBin3->Draw();
    drawLatexAdd(perfLabel,0.95,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s: %1.1f< #eta< %1.1f", calo.Data(), nominalEtaRegion[calDir][0], nominalEtaRegion[calDir][1]),0.95,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("#it{E} = %1.1f GeV", partE[exampleBin[iSet]]),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cExampleBin->Print(Form("%s/%s_%s_ExampleBinHighest.%s", outputDir.Data(), calo.Data(), readNames[pid].Data(), suffix.Data()));    
  }

  
  

    Double_t arrayBoundariesEResoX1_4[4];
    Double_t arrayBoundariesEResoY1_4[3];
    Double_t relativeMarginsEResoX[3];
    Double_t relativeMarginsEResoY[3];
    textSizeLabelsPixel             = 750*0.05;
    ReturnCorrectValuesForCanvasScaling(2250,750, 3, 1,0.04, 0.01, 0.01,0.095,arrayBoundariesEResoX1_4,arrayBoundariesEResoY1_4,relativeMarginsEResoX,relativeMarginsEResoY);

    TCanvas* canvasEResoPaperPlot     = new TCanvas("canvasEResoPaperPlot","",0,0,2250,750);  // gives the page size
    DrawGammaCanvasSettings( canvasEResoPaperPlot,  0.05, 0.01, 0.01,0.095);
    // canvasEResoPaperPlot->SetLogy(1);

    TPad* padEResoEEMC               = new TPad("padEResoEEMC", "", arrayBoundariesEResoX1_4[0], arrayBoundariesEResoY1_4[2], arrayBoundariesEResoX1_4[1], arrayBoundariesEResoY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEResoEEMC, relativeMarginsEResoX[0], relativeMarginsEResoX[1], relativeMarginsEResoY[0], relativeMarginsEResoY[2]);
    padEResoEEMC->Draw();

    TPad* padEResoBEMC                = new TPad("padEResoBEMC", "", arrayBoundariesEResoX1_4[1], arrayBoundariesEResoY1_4[2], arrayBoundariesEResoX1_4[2], arrayBoundariesEResoY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEResoBEMC, relativeMarginsEResoX[1], relativeMarginsEResoX[1], relativeMarginsEResoY[0], relativeMarginsEResoY[2]);
    padEResoBEMC->Draw();

    TPad* padEResoFEMC                = new TPad("padEResoFEMC", "", arrayBoundariesEResoX1_4[2], arrayBoundariesEResoY1_4[2], arrayBoundariesEResoX1_4[3], arrayBoundariesEResoY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEResoFEMC, relativeMarginsEResoX[1], relativeMarginsEResoX[2], relativeMarginsEResoY[0], relativeMarginsEResoY[2]);
    padEResoFEMC->Draw();

    // padEResoEEMC->cd();
    // padEResoEEMC->SetLogy();

    Double_t margin                 = relativeMarginsEResoX[0]*2.7*1350;
    double textsizeLabelsEReso    = 0;
    double textsizeFacEReso       = 0;
    if (padEResoEEMC->XtoPixel(padEResoEEMC->GetX2()) < padEResoEEMC->YtoPixel(padEResoEEMC->GetY1())){
        textsizeLabelsEReso         = (Double_t)textSizeLabelsPixel/padEResoEEMC->XtoPixel(padEResoEEMC->GetX2()) ;
        textsizeFacEReso            = (Double_t)1./padEResoEEMC->XtoPixel(padEResoEEMC->GetX2()) ;
    } else {
        textsizeLabelsEReso         = (Double_t)textSizeLabelsPixel/padEResoEEMC->YtoPixel(padEResoEEMC->GetY1());
        textsizeFacEReso            = (Double_t)1./padEResoEEMC->YtoPixel(padEResoEEMC->GetY1());
    }
    TH2F* dummyHistoEReso   = new TH2F("dummyHistoEReso","dummyHistoEReso",1000,0.001, 1.25,1000,0.0, 0.229);
    SetStyleHistoTH2ForGraphs(dummyHistoEReso, "#it{E}^{rec}/#it{E}^{MC}","probability", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,0.85*textSizeSinglePad, 0.9,1.4);
    dummyHistoEReso->GetXaxis()->SetNoExponent();
    // dummyHistoEReso->GetYaxis()->SetRangeUser(0.0,1.35);
    // dummyHistoEReso->GetXaxis()->SetRangeUser(0.15,50);
    dummyHistoEReso->GetYaxis()->SetNdivisions(505,kTRUE);
    dummyHistoEReso->GetXaxis()->SetMoreLogLabels(kTRUE);
    // histo2DPi0EResoDummy->GetYaxis()->SetTitle("norm. counts");
    // histo2DPi0EResoDummy->GetYaxis()->SetTitleOffset(1.0);
    // histo2DPi0EResoDummy->GetXaxis()->SetTickLength(0.025);
    // histo2DPi0EResoDummy->GetXaxis()->SetRangeUser(0,1.99);
    // histo2DPi0EResoDummy->GetYaxis()->SetRangeUser(0.1,2.9);
    TLegend* legendEffiE[10] = {nullptr};
    TH1D* histResoEbinshighestFitsExample[maxNSets][nPID] = {NULL};
    for (Int_t iSet = 0; iSet < 3; iSet++){
      if(iSet==1){
        padEResoBEMC->cd();
        // padEResoBEMC->SetLogx();
      } else if(iSet==2){
        padEResoFEMC->cd();
        // padEResoFEMC->SetLogx();
      } else {
        padEResoEEMC->cd();
        // padEResoEEMC->SetLogx();
      }

      dummyHistoEReso->DrawCopy();
      if(iSet==2)drawLatexAdd( Form("%s, #it{E} = %1.1f GeV", labels[iSet].Data(), partE[exampleBin[iSet]]),0.93,0.91,textSizeLabelsRel,kFALSE,kFALSE,true);
      else drawLatexAdd( Form("%s, #it{E} = %1.1f GeV", labels[iSet].Data(), partE[exampleBin[iSet]]),0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,true);

      if(iSet==0)drawLatexAdd( "a)",0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,true);
      else if(iSet==1)drawLatexAdd( "b)",0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,true);
      else drawLatexAdd( "c)",0.93,0.15,textSizeLabelsRel,kFALSE,kFALSE,true);
          // drawLatexAdd(Form("#it{E} = %1.1f GeV", partE[exampleBin[iSet]]),0.17,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);


      if(iSet==0){
        double startleg = 0.95-(nLinesCol)*textSizeLabelsRel*1.1;
        drawLatexAdd(perfLabel,0.15,startleg,textSizeLabelsRel,kFALSE,kFALSE,false);
        drawLatexAdd(Form("single particles"),0.15, startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,false);
        // drawLatexAdd(Form("#varepsilon_{TM} = N_{clus}^{matched} / N_{tracks}^{in acc.}"),0.95, startleg-(1+nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,true);
      }
      DrawGammaLines(0.15,50, 1., 1., 1, kGray+2, 7);
      TLegend* legendExamBin  = GetAndSetLegend2(0.97-0.1*nActivePartAll[iSet], 0.84, 0.97, 0.84+(textSizeLabelsRel),textSizeLabelsRel, nActivePartAll[iSet], "", 42, 0.35);
      legendExamBin->SetMargin(0.5);
      for (Int_t pid = 0; pid <nPID-1; pid++){
        if (!have_exampleBinF[iSet][pid]){
          cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << " for ExampleBinHighest plot"  << endl;
          continue;
        }

        DrawGammaSetMarker(h_exampleBinFH_reso[iSet][pid], markerStylePart[pid], markerSizePart[pid]*0.7, colorPart[pid], colorPart[pid]);
        h_exampleBinFH_reso[iSet][pid]->Draw("p,e,same");
        legendExamBin->AddEntry(h_exampleBinFH_reso[iSet][pid], labelPart[pid], "p");
        
        if (histResoEbinsFitsExampleHighest[iSet][pid] ){
          DrawGammaSetMarker(histResoEbinsFitsExampleHighest[iSet][pid], 0, 0, colorPart[pid]+1, colorPart[pid]+1, lineStylePart[pid], 3);
          histResoEbinsFitsExampleHighest[iSet][pid]->Draw("lc,same");
        }
      }
      if(iSet==0)legendExamBin->Draw();
      dummyHistoEReso->Draw("same,axis");
    }
    canvasEResoPaperPlot->SaveAs(Form("%s/EReso_PaperPlot_ECals.%s",outputDir.Data(), suffix.Data()));


    Double_t arrayBoundariesEffEtaHCalX1_4[4];
    Double_t arrayBoundariesEffEtaHCalY1_4[3];
    Double_t relativeMarginsEffEtaHCalX[3];
    Double_t relativeMarginsEffEtaHCalY[3];
    textSizeLabelsPixel             = 750*0.05;
    ReturnCorrectValuesForCanvasScaling(1500,750, 2, 1,0.065, 0.01, 0.01,0.095,arrayBoundariesEffEtaHCalX1_4,arrayBoundariesEffEtaHCalY1_4,relativeMarginsEffEtaHCalX,relativeMarginsEffEtaHCalY);

    TCanvas* canvasEffEtaHCalPaperPlot     = new TCanvas("canvasEffEtaHCalPaperPlot","",0,0,1500,750);  // gives the page size
    DrawGammaCanvasSettings( canvasEffEtaHCalPaperPlot,  0.05, 0.01, 0.01,0.095);
    // canvasEffEtaHCalPaperPlot->SetLogy(1);

    TPad* padEffEtaHCalOHCAL               = new TPad("padEffEtaHCalOHCAL", "", arrayBoundariesEffEtaHCalX1_4[0], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalOHCAL, relativeMarginsEffEtaHCalX[0], relativeMarginsEffEtaHCalX[1], relativeMarginsEffEtaHCalY[0], relativeMarginsEffEtaHCalY[2]);
    padEffEtaHCalOHCAL->Draw();

    TPad* padEffEtaHCalLFHCAL                = new TPad("padEffEtaHCalLFHCAL", "", arrayBoundariesEffEtaHCalX1_4[1], arrayBoundariesEffEtaHCalY1_4[2], arrayBoundariesEffEtaHCalX1_4[2], arrayBoundariesEffEtaHCalY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padEffEtaHCalLFHCAL, relativeMarginsEffEtaHCalX[1], relativeMarginsEffEtaHCalX[2], relativeMarginsEffEtaHCalY[0], relativeMarginsEffEtaHCalY[2]);
    padEffEtaHCalLFHCAL->Draw();


    double textsizeLabelsEffEtaHCal    = 0;
    double textsizeFacEffEtaHCal       = 0;
    if (padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) < padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1())){
        textsizeLabelsEffEtaHCal         = (Double_t)textSizeLabelsPixel/padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) ;
        textsizeFacEffEtaHCal            = (Double_t)1./padEffEtaHCalOHCAL->XtoPixel(padEffEtaHCalOHCAL->GetX2()) ;
    } else {
        textsizeLabelsEffEtaHCal         = (Double_t)textSizeLabelsPixel/padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1());
        textsizeFacEffEtaHCal            = (Double_t)1./padEffEtaHCalOHCAL->YtoPixel(padEffEtaHCalOHCAL->GetY1());
    }


    dummyHistoEReso   = new TH2F("dummyHistoEReso","dummyHistoEReso",1000,0.001, 1.75,1000,0.0, 0.109);
    SetStyleHistoTH2ForGraphs(dummyHistoEReso, "#it{E}^{rec}/#it{E}^{MC}","probability", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,0.85*textSizeSinglePad, 0.9,1.4);
    dummyHistoEReso->GetXaxis()->SetNoExponent();
    dummyHistoEReso->GetYaxis()->SetNoExponent();
    // dummyHistoEReso->GetYaxis()->SetRangeUser(0.0,1.35);
    // dummyHistoEReso->GetXaxis()->SetRangeUser(0.15,50);
    dummyHistoEReso->GetYaxis()->SetNdivisions(505,kTRUE);
    dummyHistoEReso->GetXaxis()->SetMoreLogLabels(kTRUE);
    for (Int_t iSet = 3; iSet < 5; iSet++){
      if(iSet==3){
        padEffEtaHCalOHCAL->cd();
        // padEffEtaHCalOHCAL->SetLogx();
      } else if(iSet==4){
        padEffEtaHCalLFHCAL->cd();
        // padEffEtaHCalLFHCAL->SetLogx();
      }

      dummyHistoEReso->DrawCopy();
      if(iSet==3)drawLatexAdd( Form("%s, #it{E} = %1.1f GeV", labels[iSet].Data(), partE[exampleBin[iSet]]),0.96,0.91,textSizeLabelsRel,kFALSE,kFALSE,true);
      else drawLatexAdd(Form("%s, #it{E} = %1.1f GeV", labels[iSet].Data(), partE[exampleBin[iSet]]),0.95,0.91,textSizeLabelsRel,kFALSE,kFALSE,true);
      if(iSet==3)drawLatexAdd( "d)",0.95,0.17,textSizeLabelsRel,kFALSE,kFALSE,true);
      else drawLatexAdd("e)",0.93,0.15,textSizeLabelsRel,kFALSE,kFALSE,true);

      if(iSet==3){
        double startleg = 0.95-(nLinesCol)*textSizeLabelsRel*1.1;
        drawLatexAdd(perfLabel,0.17,startleg,textSizeLabelsRel,kFALSE,kFALSE,false);
        drawLatexAdd(Form("single particles"),0.17, startleg-(nLinesCol)*textSizeLabelsRel*1.1,textSizeLabelsRel,kFALSE,kFALSE,false);
      }
      DrawGammaLines(0.15,99, 1., 1., 1, kGray+2, 7);
      TLegend* legendExamBin  = GetAndSetLegend2(0.97-0.1*nActivePartAll[iSet], 0.84, 0.97, 0.84+(textSizeLabelsRel),textSizeLabelsRel, nActivePartAll[iSet], "", 42, 0.35);
      legendExamBin->SetMargin(0.5);
      for (Int_t pid = 0; pid <nPID-1; pid++){
        if (!have_exampleBinF[iSet][pid]){
          cout << "can not find exampleBin for " << labels[iSet].Data() << " for particle " << readNames[pid] << " for ExampleBinHighest plot"  << endl;
          continue;
        }

        DrawGammaSetMarker(h_exampleBinFH_reso[iSet][pid], markerStylePart[pid], markerSizePart[pid]*0.7, colorPart[pid], colorPart[pid]);
        h_exampleBinFH_reso[iSet][pid]->Draw("p,e,same");
        legendExamBin->AddEntry(h_exampleBinFH_reso[iSet][pid], labelPart[pid], "p");

        if (histResoEbinsFitsExampleHighest[iSet][pid] ){
          DrawGammaSetMarker(histResoEbinsFitsExampleHighest[iSet][pid], 0, 0, colorPart[pid]+1, colorPart[pid]+1, lineStylePart[pid], 3);
          histResoEbinsFitsExampleHighest[iSet][pid]->Draw("lc,same");
        }
      }
      if(iSet==3)legendExamBin->Draw();
      dummyHistoEReso->Draw("same,axis");
    }
      canvasEffEtaHCalPaperPlot->SaveAs(Form("%s/EReso_PaperPlot_HCals.%s",outputDir.Data(), suffix.Data()));

  
  
}
