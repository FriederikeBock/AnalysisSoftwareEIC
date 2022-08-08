#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void positionResolutionCalorimeters_comparePlot(
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
  Style_t markerStyleSet2[8]  = {24, 46, 27,  25, 30, 42, 46, 20};
  Style_t markerStyleSet[8]   = {20, 47, 33,  21, 29, 42, 46, 20};
  Size_t markerSizeSet[8]     = {1.8*1.5, 2.0*1.5, 2.4*1.5, 1.5*1.5, 1.8*1.5, 1.8, 1.5, 1.5 };
  
  
  const Int_t maxNSets        = 10;
  
  TH1D* h_EtaReso_E[maxNSets][nPID]     = {{NULL}};
  TGraphErrors* graph_EtaReso_E[maxNSets][nPID]     = {{NULL}};
  TH1D* h_EtaMean_E[maxNSets][nPID]     = {{NULL}};
  TH1D* h_PhiReso_E[maxNSets][nPID]     = {{NULL}};
  TGraphErrors* graph_PhiReso_E[maxNSets][nPID]     = {{NULL}};
  TH1D* h_PhiMean_E[maxNSets][nPID]     = {{NULL}};

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
      colorSet[nSets] = getCaloColor(labels[nSets]);
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
      h_EtaMean_E[iSet][pid]             = (TH1D*)inputFiles[iSet]->Get(Form("h_%sreso_EtaMean_E", readNames[pid].Data()));
      if(h_EtaMean_E[iSet][pid]) currPIDavail[pid] = true;
      else cout << Form("h_%sreso_EtaMean_E", readNames[pid].Data()) << " not found" << endl;
      h_EtaReso_E[iSet][pid]             = (TH1D*)inputFiles[iSet]->Get(Form("h_%sreso_EtaSigma_E", readNames[pid].Data()));
      if(!h_EtaReso_E[iSet][pid] ) cout << Form("h_%sreso_EtaSigma_E", readNames[pid].Data()) << " not found" << endl;
      else {
        graph_EtaReso_E[iSet][pid]     = new TGraphErrors(h_EtaReso_E[iSet][pid]);
        RemoveZerosFromGraph(graph_EtaReso_E[iSet][pid]); 
      }
      h_PhiMean_E[iSet][pid]             = (TH1D*)inputFiles[iSet]->Get(Form("h_%sreso_PhiMean_E", readNames[pid].Data()));
      if(!h_PhiMean_E[iSet][pid]) cout << Form("h_%sreso_PhiMean_E", readNames[pid].Data()) << " not found" << endl;
      h_PhiReso_E[iSet][pid]             = (TH1D*)inputFiles[iSet]->Get(Form("h_%sreso_PhiSigma_E", readNames[pid].Data()));
      if(!h_PhiReso_E[iSet][pid] ) {
        cout << Form("h_%sreso_PhiSigma_E", readNames[pid].Data()) << " not found" << endl;
      } else {
        graph_PhiReso_E[iSet][pid]     = new TGraphErrors(h_PhiReso_E[iSet][pid]);
        RemoveZerosFromGraph(graph_PhiReso_E[iSet][pid]); 
      }
    }
  }

  //*******************************************************************
  // Define Canvas and pads for vertical 2 panel plots
  // -----
  // |   |
  // -----
  // |   | 
  // -----
  //*******************************************************************
  Double_t arrayBoundariesX1_4[2];
  Double_t arrayBoundariesY1_4[3];
  Double_t relativeMarginsX[3];
  Double_t relativeMarginsY[3];
  textSizeLabelsPixel             = 50;
  ReturnCorrectValuesForCanvasScaling(1350,1450, 1, 2,0.105, 0.006, 0.005,0.075,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

  TCanvas* canvasCombinedReso     = new TCanvas("canvasCombinedReso","",0,0,1350,1450);  // gives the page size
  DrawGammaCanvasSettings( canvasCombinedReso,  0.13, 0.02, 0.03, 0.06);

  TPad* padEtaReso               = new TPad("padEtaReso", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
  DrawGammaPadSettings( padEtaReso, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
  padEtaReso->Draw();
  Double_t margin                 = relativeMarginsX[0]*2.7*1350;
  Double_t textsizeLabelsECals    = 0;
  if (padEtaReso->XtoPixel(padEtaReso->GetX2()) < padEtaReso->YtoPixel(padEtaReso->GetY1())){
      textsizeLabelsECals         = (Double_t)textSizeLabelsPixel/padEtaReso->XtoPixel(padEtaReso->GetX2()) ;
  } else {
      textsizeLabelsECals         = (Double_t)textSizeLabelsPixel/padEtaReso->YtoPixel(padEtaReso->GetY1());
  }
  

  TPad* padPhiReso                = new TPad("padPhiReso", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
  DrawGammaPadSettings( padPhiReso, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
  padPhiReso->Draw();
  Double_t textsizeLabelsHCal         = 0;
  if (padPhiReso->XtoPixel(padPhiReso->GetX2()) <padPhiReso->YtoPixel(padPhiReso->GetY1()) ){
      textsizeLabelsHCal              = (Double_t)textSizeLabelsPixel/padPhiReso->XtoPixel(padPhiReso->GetX2()) ;
  } else {
      textsizeLabelsHCal              = (Double_t)textSizeLabelsPixel/padPhiReso->YtoPixel(padPhiReso->GetY1());
  }

  //*******************************************************************
  // plotting Eta reso in corresponding panels
  //*******************************************************************
  // EtaReso
  //*******************************************************************
  padEtaReso->cd();
  padEtaReso->SetLogx();

  TLegend* legendPtResMLeft2  = GetAndSetLegend2(0.5, 0.93-(3*textsizeLabelsHCal), 0.95, 0.93,textSizeLabelsPixel, 2, "", 43, 0.25);
  
  TH2F * histo2DECalEtaReso    = new TH2F("histo2DECalEtaReso","histo2DECalEtaReso",1000,0.21, 55,1000,0.00001,0.125);
  SetStyleHistoTH2ForGraphs(histo2DECalEtaReso, "#it{E} (GeV)","#sigma(#eta^{rec} - #eta^{MC})", 0.85*textsizeLabelsECals, textsizeLabelsECals,
                            0.85*textsizeLabelsECals, textsizeLabelsECals, 0.8,0.7,512,505);//#it{p}_{T} (GeV/#it{c})
  histo2DECalEtaReso->GetYaxis()->SetMoreLogLabels(kTRUE);
  // histo2DECalEtaReso->GetYaxis()->SetNdivisions(505);
  histo2DECalEtaReso->GetYaxis()->SetNoExponent(kTRUE);
  histo2DECalEtaReso->GetXaxis()->SetTickLength(0.05);
  histo2DECalEtaReso->GetYaxis()->SetTickLength(0.023);
  histo2DECalEtaReso->SetMaximum(20);
  histo2DECalEtaReso->DrawCopy();

  for(Int_t iSet=0; iSet<3;iSet++){
    DrawGammaSetMarkerTGraph(graph_EtaReso_E[iSet][1], getCaloMarker(labels[iSet],false), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
    graph_EtaReso_E[iSet][1]->Draw("same,pe");
  } 
  for(Int_t iSet=3; iSet<nSets;iSet++){
    DrawGammaSetMarkerTGraph(graph_EtaReso_E[iSet][2], getCaloMarker(labels[iSet],false), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
    graph_EtaReso_E[iSet][2]->Draw("same,p");
  }
  
  legendPtResMLeft2->AddEntry(graph_EtaReso_E[0][1],Form("%s (e)", labels[0].Data()),"p");
  legendPtResMLeft2->AddEntry(graph_EtaReso_E[3][2],Form("%s (#pi)", labels[3].Data()),"p");
  legendPtResMLeft2->AddEntry(graph_EtaReso_E[1][1],Form("%s (e)", labels[1].Data()),"p");
  legendPtResMLeft2->AddEntry(graph_EtaReso_E[4][2],Form("%s (#pi)", labels[4].Data()),"p");
  legendPtResMLeft2->AddEntry(graph_EtaReso_E[2][1],Form("%s (e)", labels[2].Data()),"p");
  legendPtResMLeft2->Draw();
  
  drawLatexAdd("#it{#bf{ECCE}} simulation",0.135,0.88,textsizeLabelsECals,kFALSE,kFALSE,false);
  drawLatexAdd(Form("single particles"),0.135,0.81,textsizeLabelsECals,kFALSE,kFALSE,false);
  drawLatexAdd(Form("a)"),0.135,0.06,textsizeLabelsECals,kFALSE,kFALSE,false);
  histo2DECalEtaReso->Draw("same,axis"); 
  
  //*******************************************************************
  // HCals
  //*******************************************************************
  padPhiReso->cd();
  padPhiReso->SetLogx();

  TH2F * histo2DHCalEtaReso            = new TH2F("histo2DHCalEtaReso","histo2DHCalEtaReso",1000,0.21, 55,1000,0,0.21);//, 100.1, 160.9);//125.1, 155.9);
  SetStyleHistoTH2ForGraphs(histo2DHCalEtaReso, "#it{E} (GeV)","#sigma(#varphi^{rec} - #varphi^{MC})", 0.85*textsizeLabelsHCal, textsizeLabelsHCal, 0.85*textsizeLabelsHCal,
                            textsizeLabelsHCal, 0.9, 0.8,512,505);
  histo2DHCalEtaReso->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DHCalEtaReso->GetYaxis()->SetNoExponent(kTRUE);
  histo2DHCalEtaReso->GetXaxis()->SetTickLength(0.047);
  histo2DHCalEtaReso->GetYaxis()->SetTickLength(0.025);
  histo2DHCalEtaReso->GetXaxis()->SetNoExponent();
  histo2DHCalEtaReso->SetMaximum(20);
  histo2DHCalEtaReso->DrawCopy();

  for(Int_t iSet=0; iSet<3;iSet++){
    DrawGammaSetMarkerTGraph(graph_PhiReso_E[iSet][1], getCaloMarker(labels[iSet],false), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
    graph_PhiReso_E[iSet][1]->Draw("same,pe");
  } 
  for(Int_t iSet=3; iSet<nSets;iSet++){
    DrawGammaSetMarkerTGraph(graph_PhiReso_E[iSet][2], getCaloMarker(labels[iSet],false), markerSizeSet[iSet], getCaloColor(labels[iSet],false), getCaloColor(labels[iSet],false));
    graph_PhiReso_E[iSet][2]->Draw("same,pe");
  }
  drawLatexAdd(Form("b)"),0.135,0.92,textsizeLabelsHCal,kFALSE,kFALSE,false);
  histo2DHCalEtaReso->Draw("same,axis"); 
  
  canvasCombinedReso->Update();
  canvasCombinedReso->Print(Form("%s/PositionResolution_Paper.%s",outputDir.Data(),suffix.Data()));

  
}
