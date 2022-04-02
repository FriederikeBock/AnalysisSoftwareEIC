#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void energyResolutionCalorimeters_comparePlot(
                            TString configInputFiles  = "file.txt",
                            TString addName           = "",
                            TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  SetPlotStyle();
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

  TH1D* h_particle_reso_E_Eta[maxNSets][nPID][nEta]     = {{{NULL}}};
  TH1D* h_particle_resoFWHM_E_Eta[maxNSets][nPID][nEta] = {{{NULL}}};
  TGraphErrors* graph_particle_reso_E_Eta[maxNSets][nPID][nEta]     = {{{NULL}}};
  TGraphErrors* graph_particle_resoFWHM_E_Eta[maxNSets][nPID][nEta] = {{{NULL}}};
  Bool_t haveEtaSlice[maxNSets][nPID][nEta]             = {{{kFALSE}}};
  
  
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
      h_particle_reso_1oEhighest[iSet][pid]             = (TH1D*)inputFiles[iSet]->Get(Form("h_%sreso_1oEhighest", readNames[pid].Data()));
      if(h_particle_reso_1oEhighest[iSet][pid]) currPIDavail[pid] = true;
      else cout << Form("h_%sreso_1oEhighest", readNames[pid].Data()) << " not found" << endl;
      h_particle_resof_1oEhighest[iSet][pid]            = (TH1D*)inputFiles[iSet]->Get(Form("h_%sresof_1oEhighest", readNames[pid].Data()));
      if(!h_particle_resof_1oEhighest[iSet][pid]) cout << Form("h_%sresof_1oEhighest", readNames[pid].Data()) << " not found" << endl;
      fit_resohighest_particle_1oE[iSet][pid]           = (TF1*)inputFiles[iSet]->Get(Form("fit_resohighest_%s1oE", readNames[pid].Data()));
      if(!fit_resohighest_particle_1oE[iSet][pid]) cout << Form("fit_resohighest_%s1oE", readNames[pid].Data()) << " not found" << endl;
      fit_resofhighest_particle_1oE[iSet][pid]          = (TF1*)inputFiles[iSet]->Get(Form("fit_resofhighest_%s1oE", readNames[pid].Data()));
      if(!fit_resofhighest_particle_1oE[iSet][pid]) cout << Form("fit_resofhighest_%s1oE", readNames[pid].Data()) << " not found" << endl;
      for (Int_t etbin = minEtaBinCaloDis[dirCal[iSet]]; etbin < maxEtaBinCaloDis[dirCal[iSet]]; etbin++){
        h_particle_reso_E_Eta[iSet][pid][etbin]         = (TH1D*)inputFiles[iSet]->Get(Form("h_%sreso_E_Eta_%d", readNames[pid].Data(), etbin));
        h_particle_resoFWHM_E_Eta[iSet][pid][etbin]     = (TH1D*)inputFiles[iSet]->Get(Form("h_%sresof_E_Eta_%d", readNames[pid].Data(), etbin));
        if (h_particle_reso_E_Eta[iSet][pid][etbin]){
          haveEtaSlice[iSet][pid][etbin] = kTRUE;
          graph_particle_reso_E_Eta[iSet][pid][etbin]     = new TGraphErrors(h_particle_reso_E_Eta[iSet][pid][etbin]);
          RemoveZerosFromGraph(graph_particle_reso_E_Eta[iSet][pid][etbin]);
          cout << partEtaCalo[etbin] << "\t" << partEtaCalo[etbin+1] << "\t" << labelPart[pid] << endl;
          cout << "reso" << endl;
          graph_particle_reso_E_Eta[iSet][pid][etbin]->Print();
          graph_particle_resoFWHM_E_Eta[iSet][pid][etbin] = new TGraphErrors(h_particle_resoFWHM_E_Eta[iSet][pid][etbin]);
          RemoveZerosFromGraph(graph_particle_resoFWHM_E_Eta[iSet][pid][etbin]);
          cout << "reso FWHM" << endl;
          graph_particle_resoFWHM_E_Eta[iSet][pid][etbin]->Print();
        }
      }
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
                               
  // *****************************************************************************************
  // Create correct bands for Yellow Report requirements
  // *****************************************************************************************
  TGraphErrors* graphReqEM[3];
  TGraphErrors* graphReq1oEEM[3];
  TGraphErrors* graphReqHad[3];
  TGraphErrors* graphReq1oEHad[3];
  for(int idir=0;idir<3;idir++){
    graphReqEM[idir]              = new TGraphErrors(nEne); 
    graphReq1oEEM[idir]           = new TGraphErrors(nEne); 
    graphReqHad[idir]             = new TGraphErrors(nEne); 
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
    graphReqEM[idir]->Sort();
    graphReq1oEEM[idir]->Sort();
    graphReqHad[idir]->Sort();
    graphReq1oEHad[idir]->Sort();
    while (graphReq1oEHad[idir]->GetX()[graphReq1oEHad[idir]->GetN()-1] > 1 )graphReq1oEHad[idir]->RemovePoint(graphReq1oEHad[idir]->GetN()-1);

    DrawGammaSetMarkerTGraphErr(  graphReqEM[idir], 0, 0, colorSet_light2[idir], colorSet_light2[idir], 2, kTRUE, colorSet_light2[idir], kFALSE);
    DrawGammaSetMarkerTGraphErr(  graphReq1oEEM[idir], 0, 0, colorSet_light2[idir], colorSet_light2[idir], 2, kTRUE, colorSet_light2[idir], kFALSE);
    DrawGammaSetMarkerTGraphErr(  graphReqHad[idir], 0, 0, colorSet_light2[idir], colorSet_light2[idir], 2, kTRUE, colorSet_light2[idir], kFALSE);
    DrawGammaSetMarkerTGraphErr(  graphReq1oEHad[idir], 0, 0, colorSet_light2[idir], colorSet_light2[idir], 2, kTRUE, colorSet_light2[idir], kFALSE);
  }

  // *****************************************************************************************
  // plot all resolutions in one plot with YR requirements
  // *****************************************************************************************
  TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
  TH2F* histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 1.5,1000,0, 200);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  TLegend* legendPtResM  = GetAndSetLegend2(0.14,addName.Contains("HCal") ? 0.95-((2*nSets)*textSizeLabelsRel) :  0.95-(2*nSets*textSizeLabelsRel), 0.55, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
  
  for (Int_t pid = 1; pid < nPID; pid++){
    cout << readNames[pid].Data() << endl;
    if(!currPIDavail[pid]) continue;
    legendPtResM->Clear();
    histResoDummy->GetYaxis()->SetRangeUser(0,21);
    if(addName.Contains("HCal"))histResoDummy->GetYaxis()->SetRangeUser(0,139);
    if(addName.Contains("HCal"))histResoDummy->GetXaxis()->SetRangeUser(0,1.02);
    histResoDummy->Draw();

    bool iindexplotted[3] = {false};
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if(!iindexplotted[dirCal[iSet]]){
        DrawGammaSetMarkerTGraphErr(  graphReq1oEEM[dirCal[iSet]], 0, 0, getCaloColor(labels[iSet],true), getCaloColor(labels[iSet],true), 2, kTRUE, getCaloColor(labels[iSet],true), kFALSE);
        DrawGammaSetMarkerTGraphErr(  graphReq1oEHad[dirCal[iSet]], 0, 0, getCaloColor(labels[iSet],true), getCaloColor(labels[iSet],true), 2, kTRUE, getCaloColor(labels[iSet],true), kFALSE);
        if(addName.Contains("HCal")){
          if(labels[iSet].Contains("OHCAL"))graphReq1oEHad[dirCal[iSet]]->SetFillStyle(3544);
          graphReq1oEHad[dirCal[iSet]]->Draw("same,e3");
        } else{
          if(labels[iSet].Contains("FEMC"))graphReq1oEEM[dirCal[iSet]]->SetFillStyle(3545);
          if(labels[iSet].Contains("BEMC"))graphReq1oEEM[dirCal[iSet]]->SetFillStyle(3754);
          if(labels[iSet].Contains("BEMC") && (!labels[iSet].Contains("Pb-Glas") || !labels[iSet].Contains("PbGl")))          
            graphReq1oEEM[dirCal[iSet]]->Draw("same,e3");

          // if(dirCal[iSet]==1) iindexplotted[2] = true;
        }
        iindexplotted[dirCal[iSet]] = true;
      }
      
    }
    TGraphErrors* dummyGraph[maxNSets] = {nullptr};
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      cout << "\t" << calo[iSet].Data() << endl;
      // sigma resolution
      if (!h_particle_reso_1oEhighest[iSet][pid]) continue;
      DrawGammaSetMarker(h_particle_reso_1oEhighest[iSet][pid], markerStyleSet[iSet], markerSizeSet[iSet], getCaloColor(labels[iSet]), getCaloColor(labels[iSet]));
      dummyGraph[iSet] = new TGraphErrors(h_particle_reso_1oEhighest[iSet][pid]);
      DrawGammaSetMarkerTGraph(dummyGraph[iSet], markerStyleSet[iSet], markerSizeSet[iSet], getCaloColor(labels[iSet]), getCaloColor(labels[iSet]));
      dummyGraph[iSet]->SetLineStyle(lineStylePID[iSet+1]);
      dummyGraph[iSet]->SetLineWidth(3);
      h_particle_reso_1oEhighest[iSet][pid]->Draw("same,p");
      DrawGammaSetMarkerTF1(fit_resohighest_particle_1oE[iSet][pid], lineStylePID[iSet+1], 3, getCaloColor(labels[iSet]));
      if(!labels[iSet].Contains("IHCAL"))fit_resohighest_particle_1oE[iSet][pid]->Draw("same");

      if(addName.Contains("HCal")){
        if(!labels[iSet].Contains("IHCAL"))legendPtResM->AddEntry(graphReq1oEHad[dirCal[iSet]], Form("YR Requirement %s",labels[iSet].Data()),"f");
      } else {
        if(!labels[iSet].Contains("PbGl"))
        legendPtResM->AddEntry(graphReq1oEEM[dirCal[iSet]], Form("YR Requirement %s",labels[iSet].Data()),"f");
      }
      if(labels[iSet].Contains("IHCAL"))legendPtResM->AddEntry(h_particle_reso_1oEhighest[iSet][pid], Form("%s",labels[iSet].Data()),"p");
      else legendPtResM->AddEntry(dummyGraph[iSet], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labels[iSet].Data(), fit_resohighest_particle_1oE[iSet][pid]->GetParameter(1),fit_resohighest_particle_1oE[iSet][pid]->GetParameter(0)),"pl");
    }
    drawLatexAdd(labelEnergy.Data(),0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("single %s", labelPart[pid].Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendPtResM->Draw();

    cReso2->Print(Form("%s/Resolution_%sFitted.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));


    //         _       _              _ _   _       _   _       _ _
    //        | |     | |            (_) | | |     | | | |     | (_)
    //   _ __ | | ___ | |_  __      ___| |_| |__   | |_| |__   | |_ _ __   ___  ___
    //  | '_ \| |/ _ \| __| \ \ /\ / / | __| '_ \  | __| '_ \  | | | '_ \ / _ \/ __|
    //  | |_) | | (_) | |_   \ V  V /| | |_| | | | | |_| |_) | | | | | | |  __/\__ \
    //  | .__/|_|\___/ \__|   \_/\_/ |_|\__|_| |_|  \__|_.__/  |_|_|_| |_|\___||___/
    //  | |
    //  |_|
    TLegend* legendPtResM2  = GetAndSetLegend2(0.14,addName.Contains("HCal") ? 0.95-((2*nSets +1)*textSizeLabelsRel) :  0.95-((2*nSets +2)*textSizeLabelsRel), 0.55, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);

    histResoDummy->GetYaxis()->SetRangeUser(0,21);
    if(addName.Contains("HCal"))histResoDummy->GetYaxis()->SetRangeUser(0,139);
    if(addName.Contains("HCal"))histResoDummy->GetXaxis()->SetRangeUser(0,1.02);
    histResoDummy->Draw();

    bool iindexplotted2[3] = {false};
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      if(!iindexplotted2[dirCal[iSet]]){
        DrawGammaSetMarkerTGraphErr(  graphReq1oEEM[dirCal[iSet]], 0, 0, getCaloColor(labels[iSet],true), getCaloColor(labels[iSet],true), 2, kTRUE, getCaloColor(labels[iSet],true), kFALSE);
        DrawGammaSetMarkerTGraphErr(  graphReq1oEHad[dirCal[iSet]], 0, 0, getCaloColor(labels[iSet],true), getCaloColor(labels[iSet],true), 2, kTRUE, getCaloColor(labels[iSet],true), kFALSE);
        if(addName.Contains("HCal")){
          if(labels[iSet].Contains("OHCAL"))graphReq1oEHad[dirCal[iSet]]->SetFillStyle(3544);
          graphReq1oEHad[dirCal[iSet]]->Draw("same,e3");
        } else{
          if(labels[iSet].Contains("FEMC"))graphReq1oEEM[dirCal[iSet]]->SetFillStyle(3545);
          if(labels[iSet].Contains("BEMC"))graphReq1oEEM[dirCal[iSet]]->SetFillStyle(3754);
          graphReq1oEEM[dirCal[iSet]]->Draw("same,e3");
        }
        iindexplotted2[dirCal[iSet]] = true;
      }
    }
    TGraphErrors* dummyGraph2[maxNSets] = {nullptr};
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      cout << "\t" << calo[iSet].Data() << endl;
      // sigma resolution
      DrawGammaSetMarker(h_particle_reso_1oEhighest[iSet][pid], markerStyleSet[iSet], markerSizeSet[iSet], getCaloColor(labels[iSet]), getCaloColor(labels[iSet]));
      dummyGraph2[iSet] = new TGraphErrors(h_particle_reso_1oEhighest[iSet][pid]);
      DrawGammaSetMarkerTGraph(dummyGraph2[iSet], markerStyleSet[iSet], markerSizeSet[iSet], getCaloColor(labels[iSet]), getCaloColor(labels[iSet]));
      dummyGraph2[iSet]->SetLineStyle(lineStylePID[iSet+1]);
      dummyGraph2[iSet]->SetLineWidth(3);
      h_particle_reso_1oEhighest[iSet][pid]->Draw("same,p");
      DrawGammaSetMarkerTF1(fit_resohighest_particle_1oE[iSet][pid], lineStylePID[iSet+1], 3, getCaloColor(labels[iSet]));
      if(!labels[iSet].Contains("IHCAL"))fit_resohighest_particle_1oE[iSet][pid]->Draw("same");

      if(addName.Contains("HCal")){
        if(!labels[iSet].Contains("IHCAL"))legendPtResM2->AddEntry(graphReq1oEHad[dirCal[iSet]], Form("YR Requirement %s",labels[iSet].Data()),"f");
      } else {
        legendPtResM2->AddEntry(graphReq1oEEM[dirCal[iSet]], Form("YR Requirement %s",labels[iSet].Data()),"f");
      }
      if(labels[iSet].Contains("IHCAL"))legendPtResM2->AddEntry(h_particle_reso_1oEhighest[iSet][pid], Form("%s",labels[iSet].Data()),"p");
      else legendPtResM2->AddEntry(dummyGraph2[iSet], Form("%s: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labels[iSet].Data(), fit_resohighest_particle_1oE[iSet][pid]->GetParameter(1),fit_resohighest_particle_1oE[iSet][pid]->GetParameter(0)),"pl");

      TF1* tempResoFunc = NULL;
      if(labels[iSet].Contains("OHCAL")){
        tempResoFunc = new TF1(Form("fresoTB_%s",labels[iSet].Data()),"[0]+[1]*x",0.1,0.8);
        tempResoFunc->SetParameters(14.5, 75);
        DrawGammaSetMarkerTF1(tempResoFunc, 9, 3, getCaloColor(labels[iSet])+1);
        tempResoFunc->Draw("same");
      }

      if(labels[iSet].Contains("EEMC")){
        tempResoFunc = new TF1(Form("fresoTB_%s",labels[iSet].Data()),"[0]+[1]*x",0.2,1.4);
        tempResoFunc->SetParameters(1., 2.);
        DrawGammaSetMarkerTF1(tempResoFunc, 9, 3, getCaloColor(labels[iSet])+1);
        tempResoFunc->Draw("same");

      }

      if(labels[iSet].Contains("BEMC")){
        tempResoFunc = new TF1(Form("fresoTB_%s",labels[iSet].Data()),"[0]+[1]*x",0.2,1.4);
        tempResoFunc->SetParameters(1.6, 2.5);
        DrawGammaSetMarkerTF1(tempResoFunc, 6, 3, getCaloColor(labels[iSet])+1);
        tempResoFunc->Draw("same");
      }
      if(tempResoFunc)legendPtResM2->AddEntry(tempResoFunc, Form("%s TB: #sigma/#it{E} =  %1.1f/#sqrt{#it{E}} #oplus %1.1f",labels[iSet].Data(), tempResoFunc->GetParameter(1),tempResoFunc->GetParameter(0)),"l");

    }
    drawLatexAdd(labelEnergy.Data(),0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("single %s", labelPart[pid].Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendPtResM2->Draw();

    cReso2->Print(Form("%s/Resolution_%sFitted_wTB.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));

  }
  
  TH2F * histResoDummyE                           = new TH2F("histResoDummyE","histResoDummyE",1000,0.3, 50,1000,0, 100);
  SetStyleHistoTH2ForGraphs(histResoDummyE, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummyE->GetXaxis()->SetNoExponent();
  histResoDummyE->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummyE->GetXaxis()->SetMoreLogLabels(kTRUE);  
  histResoDummyE->GetYaxis()->SetRangeUser(0,16.5);

  TH2F * histResoFDummyE                           = new TH2F("histResoFDummyE","histResoFDummyE",1000,0.3, 50,1000,0, 100);
  SetStyleHistoTH2ForGraphs(histResoFDummyE, "#it{E} (GeV)","FWHM / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoFDummyE->GetXaxis()->SetNoExponent();
  histResoFDummyE->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoFDummyE->GetXaxis()->SetMoreLogLabels(kTRUE);  
  histResoFDummyE->GetYaxis()->SetRangeUser(0,16.5);
  
  
  TLegend* legendResEta  = GetAndSetLegend2(0.52,0.86, 0.97, 0.86-8*0.75*textSizeLabelsRel,0.75*textSizeLabelsRel, 2, "", 42, 0.15);
  TH2D* reso2D[nPID];
  TH2D* reso2DFWHM[nPID];
  
  TCanvas* cReso3 = new TCanvas("cReso3","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso3, 0.095, 0.01, 0.01, 0.105);

  TCanvas* cReso4 = new TCanvas("cReso4","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso4, 0.11, 0.16, 0.01, 0.105);
  
  for (Int_t pid = 1; pid < nPID; pid++){
    cReso3->cd();
    cReso3->SetLogx();

    legendResEta->Clear();
    cout << readNames[pid].Data() << endl;
    if(!currPIDavail[pid]) continue;
    if (pid == 1 || pid == 0)
      histResoDummyE->GetYaxis()->SetRangeUser(0,21);
    else 
      histResoDummyE->GetYaxis()->SetRangeUser(0,70);
    histResoDummyE->Draw();

    reso2D[pid]       = new TH2D(Form("%sreso", readNames[pid].Data()), "", nEta, partEtaCalo, 2000, 0, 100);
    reso2DFWHM[pid]   = new TH2D(Form("%sresoFWHM", readNames[pid].Data()), "", nEta, partEtaCalo, 2000, 0, 100);
    
    for (Int_t iSet = 0; iSet < nSets; iSet++){
      for (Int_t iEta = 0; iEta < nEta; iEta++ ){
        if (!haveEtaSlice[iSet][pid][iEta]) continue;
        DrawGammaSetMarkerTGraph(graph_particle_reso_E_Eta[iSet][pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        graph_particle_reso_E_Eta[iSet][pid][iEta]->Draw("same,ep");
        for (Int_t k = 1; k < reso2D[pid]->GetNbinsY()+1; k++){
          if (reso2D[pid]->GetYaxis()->GetBinCenter(k) < graph_particle_reso_E_Eta[iSet][pid][iEta]->GetX()[0]){
            reso2D[pid]->Fill(partEtaCalo[iEta]+0.005, reso2D[pid]->GetYaxis()->GetBinCenter(k), graph_particle_reso_E_Eta[iSet][pid][iEta]->GetY()[0] );
          } else if (reso2D[pid]->GetYaxis()->GetBinCenter(k) > graph_particle_reso_E_Eta[iSet][pid][iEta]->GetX()[graph_particle_reso_E_Eta[iSet][pid][iEta]->GetN()-1]){
            reso2D[pid]->Fill(partEtaCalo[iEta]+0.005, reso2D[pid]->GetYaxis()->GetBinCenter(k), graph_particle_reso_E_Eta[iSet][pid][iEta]->GetY()[graph_particle_reso_E_Eta[iSet][pid][iEta]->GetN()-1] );
          } else {
            reso2D[pid]->Fill(partEtaCalo[iEta]+0.005, reso2D[pid]->GetYaxis()->GetBinCenter(k), graph_particle_reso_E_Eta[iSet][pid][iEta]->Eval(reso2D[pid]->GetYaxis()->GetBinCenter(k)) );
          }
          if (reso2DFWHM[pid]->GetYaxis()->GetBinCenter(k) < graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->GetX()[0]){
            reso2DFWHM[pid]->Fill(partEtaCalo[iEta]+0.005, reso2DFWHM[pid]->GetYaxis()->GetBinCenter(k), graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->GetY()[0] );
          } else if (reso2DFWHM[pid]->GetYaxis()->GetBinCenter(k) > graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->GetX()[graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->GetN()-1]){
            reso2DFWHM[pid]->Fill(partEtaCalo[iEta]+0.005, reso2DFWHM[pid]->GetYaxis()->GetBinCenter(k), graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->GetY()[graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->GetN()-1] );
          } else {
            reso2DFWHM[pid]->Fill(partEtaCalo[iEta]+0.005, reso2DFWHM[pid]->GetYaxis()->GetBinCenter(k), graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->Eval(reso2DFWHM[pid]->GetYaxis()->GetBinCenter(k)) );
          }
          
        }
        for (Int_t k = 1; k < reso2D[pid]->GetNbinsY()+1; k++){
          for (Int_t l = 1; l < reso2D[pid]->GetNbinsX()+1; l++){
            reso2D[pid]->SetBinError(l, k, 0);
            reso2DFWHM[pid]->SetBinError(l, k, 0);
          }
        }
        legendResEta->AddEntry(graph_particle_reso_E_Eta[iSet][pid][iEta], Form("%0.1f < #eta < %0.1f", partEtaCalo[iEta], partEtaCalo[iEta+1]), "p" );
      }
      
    }
    drawLatexAdd(labelEnergy.Data(),0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("single %s", labelPart[pid].Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendResEta->Draw();
    cReso3->Print(Form("%s/Resolution_%sAllBins.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));

    cReso3->cd();
    cReso3->SetLogx();

    if (pid == 1 || pid == 0)
      histResoFDummyE->GetYaxis()->SetRangeUser(0,21);
    else 
      histResoFDummyE->GetYaxis()->SetRangeUser(0,70);
    histResoFDummyE->Draw();

    for (Int_t iSet = 0; iSet < nSets; iSet++){
      for (Int_t iEta = 0; iEta < nEta; iEta++ ){
        if (!haveEtaSlice[iSet][pid][iEta]) continue;
        DrawGammaSetMarkerTGraph(graph_particle_resoFWHM_E_Eta[iSet][pid][iEta], markerStyleEta[iEta], markerSizeEta[iEta], colorEta[iEta], colorEta[iEta]);
        graph_particle_resoFWHM_E_Eta[iSet][pid][iEta]->Draw("same,ep");
      }
      
    }
    drawLatexAdd(labelEnergy.Data(),0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("single %s", labelPart[pid].Data()),0.95,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    legendResEta->Draw();
    cReso3->Print(Form("%s/ResolutionFWHM_%sAllBins.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));
    
    cReso4->cd();
    cReso4->SetLogy();
    cReso4->SetLogz();
      SetStyleHistoTH2ForGraphs(reso2D[pid], "#eta","#it{E} (GeV)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1);
      reso2D[pid]->SetZTitle("#sigma_{E}/E (%)");
      reso2D[pid]->GetYaxis()->SetRangeUser(0.3,99);
      reso2D[pid]->GetYaxis()->SetMoreLogLabels();
      reso2D[pid]->GetZaxis()->SetMoreLogLabels();
      reso2D[pid]->Draw("colz");
      drawLatexAdd(labelEnergy.Data(),0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("single %s", labelPart[pid].Data()),0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    
    cReso4->Print(Form("%s/Resolution_%s2D.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));
   
    SetStyleHistoTH2ForGraphs(reso2DFWHM[pid], "#eta","#it{E} (GeV)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1);
      reso2DFWHM[pid]->SetZTitle("FWHM_{E}/E (%)");
      reso2DFWHM[pid]->GetYaxis()->SetRangeUser(0.3,99);
      reso2DFWHM[pid]->GetYaxis()->SetMoreLogLabels();
      reso2DFWHM[pid]->GetZaxis()->SetMoreLogLabels();
      reso2DFWHM[pid]->Draw("colz");
      drawLatexAdd(labelEnergy.Data(),0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("single %s", labelPart[pid].Data()),0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    
    cReso4->Print(Form("%s/ResolutionFWHM_%s2D.%s", outputDir.Data(), readNames[pid].Data(), suffix.Data()));
  }
  
  TFile* fileOutput = new TFile(Form("%s/resolutionSummary.root", outputDir.Data()),"RECREATE");
  for (Int_t pid = 1; pid < nPID; pid++){
    if(!currPIDavail[pid]) continue;
    reso2D[pid]->Write();
    reso2DFWHM[pid]->Write();
  }
  fileOutput->Write();
  fileOutput->Close();
  
}
