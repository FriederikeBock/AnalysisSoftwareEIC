#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);


Double_t maxResPerEM[3]       = {8.5, 10.5, 17.5};
Double_t maxResPerHad[3]      = {62, 82, 72};
Double_t maxResEM[3]          = {1.125, 1.175, 1.205};
Double_t maxResHad[3]         = {1.25, 1.85, 1.65};

void energyResolutionCalorimeters_comparePlot(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString suffix            = "pdf"
                            
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput             = ReturnDateStringForOutput();


  TString detLabel = "";

  TString outputDir                 = Form("plots/%s/Compare%s",dateForOutput.Data(), addName.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  // for (Int_t iCl = 0; iCl < nClus; iCl++){
  //   gSystem->Exec("mkdir -p "+outputDir+"/"+calo+nameClus);
  // }
  TString labelEnergy   = "#it{#bf{ECCE}} simulation";
  TString collisionSystem = GetCollisionEnergy(addName);
  TString pTHard = "";
  Int_t nLinesCol = 1;
    
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;
  Bool_t debugOutput = kFALSE;

  Int_t nParticles                  = 7;
  TString readNames[7]              = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_", ""};
  TString labelPart[7]              = {"#gamma", "e^{#pm}", "#pi^{#pm}", "p/#bar{p}", "K^{#pm}", "n", "all"};

  Color_t colorSet[8]                   = {kBlack, kRed+1, kGreen+2, kGray+1, kCyan+1, kAzure+2, kBlack };
  Color_t colorSet_light[nPID]          = {kGray+1, kRed-6, kGreen-6, kCyan-6, kBlue+1, kOrange};
  Color_t colorSet_light2[nPID]         = {kGray+2, kRed-2, kGreen-2, kCyan-6, kBlue+1, kOrange};
  Style_t markerStyleSet[8]   = {20, 47, 33,  21, 29, 43, 46, 20};
  Style_t markerStyleSet2[8]  = {24, 46, 27,  25, 30, 42, 46, 20};
  Size_t markerSizeSet[8]     = {1.8, 2.0, 2.4, 1.5, 1.8, 1.8, 1.5, 1.5 };
  
  const Int_t maxNSets        = 10;
  
  TH1D* h_particle_reso_1oEhighest[maxNSets][nPID]     = {{NULL}};
  TH1D* h_particle_resof_1oEhighest[maxNSets][nPID]     = {{NULL}};
  TF1*  fit_resohighest_particle_1oE[maxNSets][nPID]   = {{NULL}};
  TF1*  fit_resofhighest_particle_1oE[maxNSets][nPID]   = {{NULL}};

  TString inputFilesNames[maxNSets];
  TFile* inputFiles[maxNSets];
  TString labels[maxNSets];
  TString calo[maxNSets];
  Int_t dirCal[maxNSets];

  bool currPIDavail[nPID] = {false};

  // read folder and name from file
  ifstream in(configInputFiles.Data());
  cout<<"Available Cuts:"<<endl;
  Int_t nSets = 0;
  
  while(!in.eof() ){
      in >> inputFilesNames[nSets] >> calo[nSets]>> labels[nSets]>> dirCal[nSets];
      std::cout << nSets << "\t"<< inputFilesNames[nSets].Data() << "\t"<< calo[nSets].Data()<< "\t"<< labels[nSets].Data() <<std::endl;
      labels[nSets].ReplaceAll("_"," ");
      // gSystem->Exec("mkdir -p "+outputDir+"/"+labels[nSets]+"MA");
      nSets++;
  }
  cout<<"=========================="<<endl;
  nSets--;
  std::cout<<"nSets: " << nSets <<std::endl;


  for (Int_t iSet = 0; iSet < nSets; iSet++){
    inputFiles[iSet]  = new TFile(inputFilesNames[iSet].Data());
    cout << calo[iSet].Data() << endl;
    for (Int_t pid = 1; pid < nParticles; pid++){
        h_particle_reso_1oEhighest[iSet][pid]      = (TH1D*)inputFiles[iSet]->Get(Form("h_%sreso_1oEhighest", readNames[pid].Data()));
        if(h_particle_reso_1oEhighest[iSet][pid]) currPIDavail[pid] = true;
        else cout << Form("h_%sreso_1oEhighest", readNames[pid].Data()) << " not found" << endl;
        h_particle_resof_1oEhighest[iSet][pid]      = (TH1D*)inputFiles[iSet]->Get(Form("h_%sresof_1oEhighest", readNames[pid].Data()));
        if(!h_particle_resof_1oEhighest[iSet][pid]) cout << Form("h_%sresof_1oEhighest", readNames[pid].Data()) << " not found" << endl;
        fit_resohighest_particle_1oE[iSet][pid]      = (TF1*)inputFiles[iSet]->Get(Form("fit_resohighest_%s1oE", readNames[pid].Data()));
        if(!fit_resohighest_particle_1oE[iSet][pid]) cout << Form("fit_resohighest_%s1oE", readNames[pid].Data()) << " not found" << endl;
        fit_resofhighest_particle_1oE[iSet][pid]      = (TF1*)inputFiles[iSet]->Get(Form("fit_resofhighest_%s1oE", readNames[pid].Data()));
        if(!fit_resofhighest_particle_1oE[iSet][pid]) cout << Form("fit_resofhighest_%s1oE", readNames[pid].Data()) << " not found" << endl;
    }
  }



  Double_t linTermECal[3][2]    = {{1,3}, {1,3}, {1,3}};
  Double_t linTermHCal[3][2]    = {{6,10}, {7,10}, {7,10}};
  Double_t sqrtETermECal[3][2]  = {{2,3}, {7,10}, {7,10}};
  Double_t sqrtETermHCal[3][2]  = {{45,50}, {85,100}, {35,50}};
  
  const int nEne = 34;
  const static Double_t partE[]   = {   0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 
                                        5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 
                                        13, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 50, 
                                        67.5, 75, 100, 150
                                    };
  TGraphErrors* graphReqEM[3];
  TGraphErrors* graphReq1oEEM[3];
  TGraphErrors* graphReqHad[3];
  TGraphErrors* graphReq1oEHad[3];
  for(int idir=0;idir<3;idir++){
    graphReqEM[idir]          = new TGraphErrors(nEne); 
    graphReq1oEEM[idir]          = new TGraphErrors(nEne); 
    graphReqHad[idir]          = new TGraphErrors(nEne); 
    graphReq1oEHad[idir]          = new TGraphErrors(nEne); 
    for (Int_t i= 0; i < nEne; i++){
  //     cout << "ebin:" << partE[i] << endl;
      Double_t value[4]             = {0};
      value[0]  = sqrtETermECal[idir][0]/TMath::Sqrt(partE[i])+linTermECal[idir][0];
      value[1]  = sqrtETermECal[idir][1]/TMath::Sqrt(partE[i])+linTermECal[idir][0];
      value[2]  = sqrtETermECal[idir][0]/TMath::Sqrt(partE[i])+linTermECal[idir][1];
      value[3]  = sqrtETermECal[idir][1]/TMath::Sqrt(partE[i])+linTermECal[idir][1];
      Double_t valueHad[4]             = {0};
      valueHad[0]  = sqrtETermHCal[idir][0]/TMath::Sqrt(partE[i])+linTermHCal[idir][0];
      valueHad[1]  = sqrtETermHCal[idir][1]/TMath::Sqrt(partE[i])+linTermHCal[idir][0];
      valueHad[2]  = sqrtETermHCal[idir][0]/TMath::Sqrt(partE[i])+linTermHCal[idir][1];
      valueHad[3]  = sqrtETermHCal[idir][1]/TMath::Sqrt(partE[i])+linTermHCal[idir][1];
  //     cout << value[0] << "\t" << value[1] << "\t"<< value[2] << "\t"<< value[3] << "\t"<< endl;
      Double_t max = -10000;
      Double_t min = 10000;
      for (Int_t k = 0; k < 4; k++){ 
          if (min > value[k]) min = value[k];
          if (max < value[k]) max = value[k];
      }
      graphReqEM[idir]->SetPoint(i, partE[i], (min+max)/2 );
      graphReqEM[idir]->SetPointError(i, 0.05, TMath::Abs(min-max)/2 );
      graphReq1oEEM[idir]->SetPoint(i, 1/TMath::Sqrt(partE[i]), (min+max)/2 );
      graphReq1oEEM[idir]->SetPointError(i, 0.005, TMath::Abs(min-max)/2 );
      
      max = -10000;
      min = 10000;
      for (Int_t k = 0; k < 4; k++){ 
          if (min > valueHad[k]) min = valueHad[k];
          if (max < valueHad[k]) max = valueHad[k];
      }
      graphReqHad[idir]->SetPoint(i, partE[i], (min+max)/2 );
      graphReqHad[idir]->SetPointError(i, 0.05, TMath::Abs(min-max)/2 );
      graphReq1oEHad[idir]->SetPoint(i, 1/TMath::Sqrt(partE[i]), (min+max)/2 );
      graphReq1oEHad[idir]->SetPointError(i, 0.005, TMath::Abs(min-max)/2 );
    }
    graphReq1oEEM[idir]->Sort();
    DrawGammaSetMarkerTGraphErr(  graphReqEM[idir], 0, 0, kBlue-10, kBlue-10, 1, kTRUE, kBlue-10, kFALSE);
    DrawGammaSetMarkerTGraphErr(  graphReq1oEEM[idir], 0, 0, kBlue-10, kBlue-10, 1, kTRUE, kBlue-10, kFALSE);
    DrawGammaSetMarkerTGraphErr(  graphReqHad[idir], 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
    DrawGammaSetMarkerTGraphErr(  graphReq1oEHad[idir], 0, 0, kGray, kGray, 1, kTRUE, kGray, kFALSE);
  }




  TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
  TH2F* histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 1.5,1000,0, 50);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  // TLegend* legendPtResM  = GetAndSetLegend2(0.14, 0.72-(nSets*textSizeLabelsRel), 0.34, 0.72,textSizeLabelsPixel, 1, "", 43, 0.15);
  TLegend* legendPtResM  = GetAndSetLegend2(0.65, 0.14, 0.95, 0.14+(nSets*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);
  

  for (Int_t pid = 1; pid < nPID; pid++){
    cout << readNames[pid].Data() << endl;
    if(!currPIDavail[pid]) continue;
    // histoDummyEffiE->Draw();
      legendPtResM->Clear();
      if(pid==1)legendPtResM  = GetAndSetLegend2(0.74, 0.15, 0.94, 0.15+(nSets*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

      histResoDummy->GetYaxis()->SetRangeUser(0,maxResPerEM[2]);
      if(addName.Contains("HCal"))histResoDummy->GetYaxis()->SetRangeUser(0,maxResPerHad[2]);
      histResoDummy->Draw();

    TLegend* legendResLabE  = GetAndSetLegend2(0.14, 0.95-0.85*textSizeLabelsRel, 0.35, 0.95,0.85*textSizeLabelsRel, 1, "", 42, 0.2);

    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if(addName.Contains("HCal"))graphReq1oEHad[dirCal[iSet]]->Draw("same,e3");
      else graphReq1oEEM[dirCal[iSet]]->Draw("same,e3");
    }
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      cout << "\t" << calo[iSet].Data() << endl;
// sigma resolution
      DrawGammaSetMarker(h_particle_reso_1oEhighest[iSet][pid], markerStyleSet[iSet], markerSizeSet[iSet], colorSet[iSet], colorSet[iSet]);
      h_particle_reso_1oEhighest[iSet][pid]->Draw("same,p");
      DrawGammaSetMarkerTF1(fit_resohighest_particle_1oE[iSet][pid], 3, 3, colorSet[iSet]);
      fit_resohighest_particle_1oE[iSet][pid]->Draw("same");

      // fwhm resolution
      DrawGammaSetMarker(h_particle_resof_1oEhighest[iSet][pid], markerStyleSet2[iSet], markerSizeSet[iSet], colorSet_light[iSet], colorSet_light[iSet]);
      h_particle_resof_1oEhighest[iSet][pid]->Draw("same,p");
      // DrawGammaSetMarkerTF1(fit_resofhighest_particle_1oE[iSet][pid], 3, 3, colorSet_light[iSet]);
      // fit_resofhighest_particle_1oE[iSet][pid]->Draw("same");
      // legendResLabE->AddEntry(fit_resohighest_particle_1oE[iSet], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[pid].Data(), fit_resohighest_particle_1oE[iSet]->GetParameter(1),fit_resohighest_particle_1oE[iSet]->GetParameter(0)),"l");
      
    }    
    // drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%s in %s", labelPart[pid].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendResLabE->Draw();

    cReso2->Print(Form("%s/Resolution_%sFitted.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));

  }




}


// TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
//     cReso->SetLogx(1);
//     DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
//     TH2F * histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0.3, 170,1000,0, maxResPerHad[dirCal]);
//     SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
//     histResoDummy->GetXaxis()->SetNoExponent();
//     histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
//     histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
//     if (i < 2) 
//       histResoDummy->GetYaxis()->SetRangeUser(0, maxResPerEM[dirCal]);
//     else 
//       histResoDummy->GetYaxis()->SetRangeUser(0, maxResPerHad[dirCal]);
//     histResoDummy->Draw();

//     TLegend* legendResLabE  = GetAndSetLegend2(0.14, 0.95-0.85*textSizeLabelsRel, 0.35, 0.95,0.85*textSizeLabelsRel, 1, "", 42, 0.2);

//     DrawGammaSetMarker(h_particle_reso_1oEhighest[iSet][pid], markerStylePart[pid], markerSizePart[pid], colorPart[pid], colorPart[pid]);
//     h_particle_reso_1oEhighest[iSet][pid]->Draw("same,p");
    
//     fit_resohighest_particle_1oE[iSet] = CreateFitE (h_particle_reso_1oEhighest[iSet][pid], Form("fit_reso_%sE", readNames[pid].Data()),resFitEM[dirCal][0], resFitEM[dirCal][1]);
//     DrawGammaSetMarkerTF1(fit_resohighest_particle_1oE[iSet], 3, 3, colorPart[pid]+1);
//     fit_resohighest_particle_1oE[iSet]->SetRange(minEPart[pid], maxEPart[pid]);
//     fit_resohighest_particle_1oE[iSet]->Draw("same");
//     if (i < 2){
//       DrawGammaSetMarker(h_particle_resof_1oEhighest[iSet][pid], markerStylePart[pid], markerSizePart[pid], colorPart[pid]-6, colorPart[pid]-6);
//       h_particle_resof_1oEhighest[iSet][pid]->Draw("same,p");
//       fit_resofhighest_particle_1oE[iSet][pid] = CreateFitE (h_particle_resof_1oEhighest[iSet][pid], Form("fit_resof_%sE", readNames[pid].Data()),resFitEM[dirCal][0], resFitEM[dirCal][1]);
//       DrawGammaSetMarkerTF1(fit_resofhighest_particle_1oE[iSet][pid], 3, 3, colorPart[pid]-6);
//       fit_resofhighest_particle_1oE[iSet][pid]->SetRange(minEPart[pid], maxEPart[pid]);
//       fit_resofhighest_particle_1oE[iSet][pid]->Draw("same");
      
//     }
    
//     drawLatexAdd(perfLabel,0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
//     drawLatexAdd(Form("%s in %s", labelPart[pid].Data() ,detLabel.Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
//     legendResLabE->AddEntry(fit_resohighest_particle_1oE[iSet], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labelPart[pid].Data(), fit_resohighest_particle_1oE[iSet]->GetParameter(1),fit_resohighest_particle_1oE[iSet]->GetParameter(0)),"l");
//     legendResLabE->Draw();

//     cReso->Print(Form("%s/Resolution_%sFitted.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));