#include <algorithm>
#include "../common/plottingheader.h"
#include "../treeAnalysis/treeProcessing.h"
#include "../treeAnalysis/clusterizer.cxx"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
#include "TColor.h"
#include <TChain.h>
#include <TVector3.h>
#include <vector>

typedef struct {
  float tower_Eta;
  float tower_Phi;
  float tower_x;
  float tower_y;
  float tower_z;
  float tower_R;
  float tower_theta;
} geomStrct;
  

void plot3DTowerstree(
    TString nameFile          = "",
    TString inFileGeometry    = "",
    TString suffix            = "pdf",
    TString collisionsSys     = "PythiaMB",
    double energy             = -1.0,
    int plotEvent             = 7,
    TString addName           = ""           
){
  int clusterizerEnum =0;
  int caloID        = 8;
  // dimensions 
  int towersx       = 105;
  int towersy       = 105;
  int towersz       = 7;
  float minECutOff  = aggE[caloID]/2.;
  float seedEC      = seedE[caloID];
  float aggEC       = aggE[caloID];
  int verbosity     = 0;
  TString labelDet  = Form("LFHCal%s", addName.Data());
  TString caloName  = "LFHCAL";
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  
  Color_t colorCluster[50] = {kRed, kGreen+2, kBlue, kOrange, kCyan+1, kMagenta+1, kViolet, kTeal, kSpring-1, kYellow,
                             kRed+3,  kGreen+3, kBlue+3, kOrange+2, kCyan+3, kMagenta+3, kViolet+2, kTeal-1, kSpring+2, kYellow+2,
                             kRed-4, kGreen-4, kBlue-4, kOrange-5, kCyan-5, kMagenta-3, kViolet-3, kTeal-6, kSpring-6, kYellow-6,
                             kRed-8, kGreen-6, kBlue-7, kOrange+7, kCyan-6, kMagenta-6, kViolet-7, kTeal+3, kSpring+9, kYellow-5,
                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  TString labelEnergy = "";
  if (collisionsSys.CompareTo("PythiaMB") == 0){
    labelEnergy   = "e-p: 18#times 275 GeV";
  } else if (collisionsSys.CompareTo("PythiaQ2_100") == 0){
    labelEnergy   = "e-p: 18#times 275 GeV, Q^{2} > 100 GeV";
  } else if (collisionsSys.CompareTo("Satre18x108") == 0){
    labelEnergy   = "e-p: 18#times 108 GeV, Satre";
  } else if (collisionsSys.CompareTo("SingleElectron") == 0){
    labelEnergy   = Form("e^{-}, %0.1f GeV", energy);
  } else if (collisionsSys.CompareTo("SingleKaon") == 0){
    labelEnergy   = Form("K^{-}, %0.1f GeV", energy);
  } else if (collisionsSys.CompareTo("SinglePion") == 0){
    labelEnergy   = Form("#pi^{-}, %0.1f GeV", energy);
  } else if (collisionsSys.CompareTo("SingleProton") == 0){
    labelEnergy   = Form("p, %0.1f GeV", energy);
  } else if (collisionsSys.CompareTo("SingleNeutron") == 0){
    labelEnergy   = Form("n, %0.1f GeV", energy);
  }
  
  TString dateForOutput                         = ReturnDateStringForOutput();

  TString outputDir                             = Form("plots/%s_%s_%s",dateForOutput.Data(),collisionsSys.Data(),caloName.Data());
  if (collisionsSys.Contains("Single"))
    outputDir = Form("%s_%.0f", outputDir.Data(), energy);
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDir+"/2D");
  gSystem->Exec("mkdir -p "+outputDir+"/2DwClusters");
  gSystem->Exec("mkdir -p "+outputDir+"/3D");
  gSystem->Exec("mkdir -p "+outputDir+"/3DwClusters");

  double textSizeSinglePad                    = 0.05;
  double textSizeLabelsPixel                  = 35;
  double textSizeLabelsRel                    = 58./1300;

  int granularity = 1;
  TH2F* h_IEtaIPhiMapEvt                   = new TH2F("h_IEtaIPhiMapEvt", "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
  TH2F* h_IEtaIPhiMapEvt_Ecut              = new TH2F("h_IEtaIPhiMapEvt_Ecut", "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
  TH2F* h_IEtaIPhiMapEvt_Hole              = new TH2F("h_IEtaIPhiMapEvt_Hole", "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
  TH2F* h_IEtaIPhiMapEvt_Cl[50]            = {NULL};
  TH2F* h_IEtaIPhiMapEvt_MCPart[50]        = {NULL};
  TH3F* h_IEtaIPhiILMapEvt                 = new TH3F("h_IEtaIPhiILMapEvt", "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5  );
  TH3F* h_IEtaIPhiILMapEvt_Hole            = new TH3F("h_IEtaIPhiILMapEvt_Hole", "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5  );
  TH3F* h_IEtaIPhiILMapEvt_Ecut            = new TH3F("h_IEtaIPhiILMapEvt_Ecut", "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
  TH3F* h_IEtaIPhiILMapEvt_Cl[50]          = {NULL};
  TH3F* h_IEtaIPhiILMapEvt_MCPart[50]      = {NULL};
  
  // read file
  // load event tree
  TChain *tt_event  = new TChain("event_tree");
  tt_event->Add(nameFile.Data());
  if( !tt_event ) { std::cout << "tree not found... returning!" << std::endl; return; }
  
  // load geometry 
  tt_geometry   = (TTree*) (new TFile(inFileGeometry.Data(),"READ"))->Get("geometry_tree");
  if( !tt_geometry ) { std:: cout << "geometry tree not found... returning! " << std:: endl; return; }

  // load all branches (see treeProcessing header)
  SetBranchAddressesTree(tt_event);
  SetBranchAddressesGeometryTree(tt_geometry);
  SetGeometryIndices();

    // load geometry of chosen calorimeter
  tt_geometry->GetEvent( calogeomindex[caloID] );

  // vector geometry structure
  geomStrct   arr_GeoStrt[towersx][towersy][towersz];  // array[iEta-1][iPhi-1]
  for(int i=0; i<(int)_calogeom_towers_N; i++){
    int tmpEta  = _calogeom_towers_iEta[i]-1;
    int tmpPhi  = _calogeom_towers_iPhi[i]-1;
    int tmpL  = _calogeom_towers_iL[i];
    h_IEtaIPhiMapEvt_Hole->SetBinContent(tmpEta+2,tmpPhi+2, 1);
    h_IEtaIPhiILMapEvt_Hole->SetBinContent(tmpL+1,tmpEta+2,tmpPhi+2, 1);
    arr_GeoStrt[tmpEta][tmpPhi][tmpL].tower_Eta   = _calogeom_towers_Eta[i];
    arr_GeoStrt[tmpEta][tmpPhi][tmpL].tower_Phi   = _calogeom_towers_Phi[i];
    arr_GeoStrt[tmpEta][tmpPhi][tmpL].tower_x     = _calogeom_towers_x[i];
    arr_GeoStrt[tmpEta][tmpPhi][tmpL].tower_y     = _calogeom_towers_y[i];
    arr_GeoStrt[tmpEta][tmpPhi][tmpL].tower_z     = _calogeom_towers_z[i];
  }

  for (int x = 1; x < h_IEtaIPhiMapEvt_Hole->GetNbinsX()+1; x++){
    for (int y = 1; y < h_IEtaIPhiMapEvt_Hole->GetNbinsY()+1; y++){
      if (h_IEtaIPhiMapEvt_Hole->GetBinContent(x,y) > 0 )
        h_IEtaIPhiMapEvt_Hole->SetBinContent(x,y, 0);
      else 
        h_IEtaIPhiMapEvt_Hole->SetBinContent(x,y, 0.1);
    }
  }

  for (int x = 1; x < h_IEtaIPhiILMapEvt_Hole->GetNbinsX()+1; x++){
    for (int y = 1; y < h_IEtaIPhiILMapEvt_Hole->GetNbinsY()+1; y++){
      for (int z = 1; z < h_IEtaIPhiILMapEvt_Hole->GetNbinsZ()+1; z++){
        if (h_IEtaIPhiILMapEvt_Hole->GetBinContent(x,y,z) > 0 )
          h_IEtaIPhiILMapEvt_Hole->SetBinContent(x,y,z, 0);
        else 
          h_IEtaIPhiILMapEvt_Hole->SetBinContent(x,y,z, 0.1);
      }
    }
  }
  
  
  float sumedE                    = 0;
  if(tt_event){
    h_IEtaIPhiMapEvt->Sumw2();
    h_IEtaIPhiMapEvt_Ecut->Sumw2();
    h_IEtaIPhiILMapEvt->Sumw2();
    h_IEtaIPhiILMapEvt_Ecut->Sumw2();
    for (int i = 0; i < 50; i++){
      h_IEtaIPhiMapEvt_Cl[i]  = new TH2F(Form("h_IEtaIPhiMapEvt_Cl%d",i), "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
      h_IEtaIPhiMapEvt_Cl[i]->Sumw2();
      h_IEtaIPhiMapEvt_MCPart[i]  = new TH2F(Form("h_IEtaIPhiMapEvt_MCPart%d",i), "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
      h_IEtaIPhiMapEvt_MCPart[i]->Sumw2();
      h_IEtaIPhiILMapEvt_Cl[i]  = new TH3F(Form("h_IEtaIPhiILMapEvt_Cl%d",i), "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
      h_IEtaIPhiILMapEvt_Cl[i]->Sumw2();
      h_IEtaIPhiILMapEvt_MCPart[i]  = new TH3F(Form("h_IEtaIPhiILMapEvt_MCPart%d",i), "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
      h_IEtaIPhiILMapEvt_MCPart[i]->Sumw2();
    }
    
    Int_t nentries = Int_t(tt_event->GetEntries());
    for (Int_t i = 0; i < nentries; ++i) {
      
      if (tt_event->LoadTree(i) < 0)
        break;
      if(i!=plotEvent) continue;
      tt_event->GetEntry(i);
      
      for(Int_t imc=0; (Int_t)imc<_nMCPart; imc++){
        TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
        _mcpart_Eta[imc]=truevec.Eta();
        _mcpart_Phi[imc]=truevec.Phi();
      }
      if(_nTowers_LFHCAL == 0 ) return;

      for(int itwr=0;(int)itwr<_nTowers_LFHCAL;itwr++){
        sumedE += _tower_LFHCAL_E[itwr];
        if (_tower_LFHCAL_E[itwr] > 1e-3 && verbosity > 2) cout << _tower_LFHCAL_iEta[itwr] << "\t" << _tower_LFHCAL_iPhi[itwr] << "\t"<< _tower_LFHCAL_iL[itwr] << "\t"<< _tower_LFHCAL_E[itwr] << endl;
        h_IEtaIPhiMapEvt->Fill(_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr],_tower_LFHCAL_E[itwr]);
        h_IEtaIPhiILMapEvt->Fill(_tower_LFHCAL_iL[itwr],_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr], _tower_LFHCAL_E[itwr]);
        if(_tower_LFHCAL_E[itwr]>0.001){
          h_IEtaIPhiMapEvt_Ecut->Fill(_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr],_tower_LFHCAL_E[itwr]);
          h_IEtaIPhiILMapEvt_Ecut->Fill(_tower_LFHCAL_iL[itwr],_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr], _tower_LFHCAL_E[itwr]);
        }  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
      }
    }
  }

  int nclusters = 0;

  // vector of towers used in the clusterizer
  std::vector<towersStrct> input_towers;
  std::vector<towersStrct> input_towersForMC;
  // vector of towers within the currently found cluster
  std::vector<towersStrct> cluster_towers;
  std::vector<clustersStrct> rec_clusters;
  std::vector<float> recclusterEnergies;
  std::vector<TString> mcPartEnergies;
  std::vector<int> clslabels;

  // fill vector with towers for clusterization above aggregation threshold
  for(int itow=0; itow<(int)_nTowers_LFHCAL; itow++){
    if(_tower_LFHCAL_E[itow]>aggEC){
      towersStrct tempstructT;
      tempstructT.tower_E       = _tower_LFHCAL_E[itow];
      tempstructT.tower_iEta    = _tower_LFHCAL_iEta[itow];
      tempstructT.tower_iPhi    = _tower_LFHCAL_iPhi[itow];
      tempstructT.tower_iL      = _tower_LFHCAL_iL[itow];
      tempstructT.tower_trueID  = _tower_LFHCAL_trueID[itow];
      input_towers.push_back(tempstructT);
      input_towersForMC.push_back(tempstructT);
      if (verbosity > 1) std::cout << tempstructT.tower_E << "\t" << tempstructT.tower_iEta << "\t" << tempstructT.tower_iPhi << "\t" << tempstructT.tower_iL << std::endl;
    }
  }

  // sort vector in descending energy order
  std::sort(input_towers.begin(), input_towers.end(), &acompare);
  std::sort(input_towersForMC.begin(), input_towersForMC.end(), &acompareTrueID);
  while (!input_towers.empty()) {
    cluster_towers.clear();
    clustersStrct tempstructC;
    clslabels.clear();

    if(input_towers.at(0).tower_E > seedEC){

      if(clusterizerEnum==kC3){
        tempstructC = findCircularCluster(3, seedEC, caloID, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kC5){
        tempstructC = findCircularCluster(5, seedEC, caloID, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==k3x3){
        tempstructC = findSquareCluster(3, seedEC, caloID, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==k5x5){
        tempstructC = findSquareCluster(5, seedEC, caloID, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kV1){  
        tempstructC = findV1Cluster(seedEC, caloID, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kV3){  
        tempstructC = findV3Cluster(seedEC, caloID, input_towers, cluster_towers, clslabels);
      } else if (clusterizerEnum==kMA){  
        tempstructC = findMACluster(seedEC, caloID, input_towers, cluster_towers, clslabels);
      } else {
        std::cout << "incorrect clusterizer selected! Enum " << clusterizerEnum << " not defined!" << std::endl;
        return;
      }
      tempstructC.cluster_towers = cluster_towers;
      cout << "filling hist for cluster: "  << nclusters << "\t E: " << tempstructC.cluster_E << "\t seed E: " << tempstructC.cluster_seed<<  endl;
      for (int tic = 0; tic < (int)cluster_towers.size(); tic++){
        int iEtaTwr = cluster_towers.at(tic).tower_iEta;
        int iPhiTwr = cluster_towers.at(tic).tower_iPhi;
        int iLTwr   = cluster_towers.at(tic).tower_iL;
        float eTwr    = cluster_towers.at(tic).tower_E;
        if (verbosity > 0) std::cout << "\t" << Form("%.3f",cluster_towers.at(tic).tower_E) << "\t eta \t" << cluster_towers.at(tic).tower_iEta << "\t phi \t" << cluster_towers.at(tic).tower_iPhi << "\t L \t" << cluster_towers.at(tic).tower_iL << std::endl;      
        h_IEtaIPhiMapEvt_Cl[nclusters]->Fill(iEtaTwr, iPhiTwr, eTwr);
        h_IEtaIPhiILMapEvt_Cl[nclusters]->Fill(iLTwr, iEtaTwr, iPhiTwr, eTwr);
      }
      rec_clusters.push_back(tempstructC);
//       recclusterEnergies.push_back(tempstructC.cluster_E);

      nclusters++;
      
    } else {
      input_towers.clear();
    }
  }
  std::sort(rec_clusters.begin(), rec_clusters.end(), &acompareCl);
  
  
    
  float maxSeedEnergy = 0;
  for (int it = 0; it < (int)rec_clusters.size(); it++){
    for (int tic = 0; tic < (int)(rec_clusters.at(it).cluster_towers).size(); tic++){
      int iEtaTwr = rec_clusters.at(it).cluster_towers.at(tic).tower_iEta;
      int iPhiTwr = rec_clusters.at(it).cluster_towers.at(tic).tower_iPhi;
      int iLTwr   = rec_clusters.at(it).cluster_towers.at(tic).tower_iL;
      float eTwr    = rec_clusters.at(it).cluster_towers.at(tic).tower_E;
      if (verbosity > 0) std::cout << "\t" << Form("%.3f",eTwr) << "\t eta \t" << iEtaTwr << "\t phi \t" << iPhiTwr << "\t L \t" << iLTwr << std::endl;      
      h_IEtaIPhiMapEvt_Cl[nclusters]->Fill(iEtaTwr, iPhiTwr, eTwr);
      h_IEtaIPhiILMapEvt_Cl[nclusters]->Fill(iLTwr, iEtaTwr, iPhiTwr, eTwr);
    }
    if (rec_clusters.at(it).cluster_seed > maxSeedEnergy)
      maxSeedEnergy   =  rec_clusters.at(it).cluster_seed;
  }
  
  if (sumedE < 0.01){
    cout << "total E calo: " << sumedE <<  " No particle hit the calorimeter! Aborting." << endl; 
    return;
  } else {
    cout << "total E calo: " << sumedE << endl; 
  }
  
  int nMCCurr   = 0;
  int nMCCurrID = input_towersForMC.at(0).tower_trueID;

  std::vector<towersStrct> cluster_towersMC;
  std::vector<clustersStrct> MC_clusters;
  clustersStrct tempstructMC;

  while (!input_towersForMC.empty() ) {
    if (nMCCurrID != input_towersForMC.at(0).tower_trueID ){
      tempstructMC.cluster_Eta    = _mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)] ;
      tempstructMC.cluster_Phi    = _mcpart_Phi[GetCorrectMCArrayEntry(nMCCurrID)] ;
      tempstructMC.cluster_trueID = nMCCurrID;
      tempstructMC.cluster_NTowers  = cluster_towersMC.size();
      tempstructMC.cluster_towers   = cluster_towersMC;
      MC_clusters.push_back(tempstructMC);
      nMCCurrID = input_towersForMC.at(0).tower_trueID;
      // reseting tempstructMC
      tempstructMC.cluster_E      = 0;
      tempstructMC.cluster_Eta    = -1000;
      tempstructMC.cluster_Phi    = -1000;
      tempstructMC.cluster_trueID = nMCCurrID;
      tempstructMC.cluster_NTowers  = 0;
      cluster_towersMC.clear();
    }
    cluster_towersMC.push_back(input_towersForMC.at(0));
    tempstructMC.cluster_E += input_towersForMC.at(0).tower_E;
    
    if (input_towersForMC.size() == 1){
      tempstructMC.cluster_Eta    = _mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)] ;
      tempstructMC.cluster_Phi    = _mcpart_Phi[GetCorrectMCArrayEntry(nMCCurrID)] ;
      tempstructMC.cluster_trueID = nMCCurrID;
      tempstructMC.cluster_NTowers  = cluster_towersMC.size();
      tempstructMC.cluster_towers   = cluster_towersMC;
      MC_clusters.push_back(tempstructMC);    
    }
    input_towersForMC.erase(input_towersForMC.begin());
  }
  std::sort(MC_clusters.begin(), MC_clusters.end(), &acompareCl);
  
  for (int nMC = 0; nMC < (int)MC_clusters.size() && nMC < 50; nMC++  ){
    int mcID    = MC_clusters.at(nMC).cluster_trueID;
    TString tempMCLabel = Form("%i: r. %.1f GeV, #eta = %1.2f", _mcpart_PDG[GetCorrectMCArrayEntry(mcID)] , MC_clusters.at(nMC).cluster_E,  MC_clusters.at(nMC).cluster_Eta  );
    mcPartEnergies.push_back(tempMCLabel);
    cout << tempMCLabel.Data() << endl;
    for (int tic = 0; tic < (int)(MC_clusters.at(nMC).cluster_towers).size(); tic++){
      int iEtaTwr = MC_clusters.at(nMC).cluster_towers.at(tic).tower_iEta;
      int iPhiTwr = MC_clusters.at(nMC).cluster_towers.at(tic).tower_iPhi;
      int iLTwr   = MC_clusters.at(nMC).cluster_towers.at(tic).tower_iL;
      float eTwr  = MC_clusters.at(nMC).cluster_towers.at(tic).tower_E;
      
      h_IEtaIPhiMapEvt_MCPart[nMC]->Fill(iEtaTwr, iPhiTwr, eTwr);
      h_IEtaIPhiILMapEvt_MCPart[nMC]->Fill(iLTwr, iEtaTwr, iPhiTwr, eTwr);
    }  
    nMCCurr++;
  }
  
  if (collisionsSys.Contains("Single"))
    labelEnergy=Form("E = %1.2f, #eta = %1.2f",_mcpart_E[0], _mcpart_Eta[0]);
  
  // 2D PLOT
  TCanvas* c2DBox = new TCanvas("c2DBox","",0,0,1000,950);
  DrawGammaCanvasSettings( c2DBox, 0.11, 0.1, 0.01, 0.1);

  // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_Ecut, "tower x index","tower y index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
    
    h_IEtaIPhiMapEvt_Ecut->Draw("colz");
    h_IEtaIPhiMapEvt_Hole->SetLineColor(kGray);
    h_IEtaIPhiMapEvt_Hole->Draw("box,same");
    
    drawLatexAdd(labelDet.Data(),0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(labelEnergy.Data(),0.88,0.15+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  c2DBox->Print(Form("%s/2D/2D_Event%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));

  // reset 2D canvas
  DrawGammaCanvasSettings( c2DBox, 0.11, 0.03, 0.01, 0.1);
  
  //************************************************************************************
  // draw 2D hist with box option
  //************************************************************************************  
    SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_Hole, "tower x index","tower y index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
    h_IEtaIPhiMapEvt_Hole->GetZaxis()->SetRangeUser(0,h_IEtaIPhiMapEvt_Ecut->GetMaximum()*1.1);
    h_IEtaIPhiMapEvt_Hole->Draw("box");
    h_IEtaIPhiMapEvt_Hole->Draw("same,axis");
    h_IEtaIPhiMapEvt_Ecut->SetLineColor(kGray+2);
    h_IEtaIPhiMapEvt_Ecut->Draw("box,same");
    
    drawLatexAdd(labelDet.Data(),0.95,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(labelEnergy.Data(),0.95,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  c2DBox->Print(Form("%s/2D/2D_Event%d_Ecut_Box.%s", outputDir.Data(),plotEvent, suffix.Data()));
  //************************************************************************************
  // draw 2D hist with box option and different clusters on top
  //************************************************************************************
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IEtaIPhiMapEvt_Cl[nCl]->SetLineColor(colorCluster[nCl]);
      h_IEtaIPhiMapEvt_Cl[nCl]->Draw("box,same");
      int cols = nCl%4;
      int rows = (nCl/4);
      drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", rec_clusters.at(nCl).cluster_E ), 0.14+0.2*cols, 0.93-(rows*0.55*textSizeLabelsRel), 0.55*textSizeLabelsRel, kFALSE, kFALSE, kFALSE, colorCluster[nCl]);
    }
    if (collisionsSys.Contains("Single")){
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.14,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
  c2DBox->Print(Form("%s/2DwClusters/2D_Event%d_Ecut_clusters.%s", outputDir.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 2D hist with box option and different MC particles on top
  //************************************************************************************  
  h_IEtaIPhiMapEvt_Hole->Draw("box");
  h_IEtaIPhiMapEvt_Ecut->SetLineColor(kGray+2);
  h_IEtaIPhiMapEvt_Ecut->Draw("same,box");
  
  for (int nMC = 0; nMC < nMCCurr; nMC++){
    if (h_IEtaIPhiMapEvt_MCPart[nMC]){
      h_IEtaIPhiMapEvt_MCPart[nMC]->SetLineColor(colorCluster[nMC]);
      h_IEtaIPhiMapEvt_MCPart[nMC]->Draw("box,same");
      int cols = nMC%4;
      int rows = (nMC/4);
      drawLatexAdd((TString)(mcPartEnergies.at(nMC)).Data() , 0.14+0.2*cols, 0.93-(rows*0.45*textSizeLabelsRel), 0.45*textSizeLabelsRel, kFALSE, kFALSE, kFALSE, colorCluster[nMC]);
    }
  }
  
  drawLatexAdd(labelDet.Data(),0.95,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(labelEnergy.Data(),0.95,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  c2DBox->Print(Form("%s/2DwClusters/2D_Event%d_Ecut_MCPart.%s", outputDir.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 3D hist with box option   
  //************************************************************************************
  TCanvas* c3DBox = new TCanvas("c3DBox","",0,0,1000,950);
  DrawGammaCanvasSettings( c3DBox, 0.12, 0.1, 0.01, 0.18);
  c3DBox->SetTheta(16);
  c3DBox->SetPhi(68);
   
    SetStyleHistoTH3ForGraphs(h_IEtaIPhiILMapEvt_Hole, "tower z index", "tower x index","tower y index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1, 1.1, 1.15, 505, 510,510);

    h_IEtaIPhiILMapEvt_Hole->SetLineColor(kGray);
    h_IEtaIPhiILMapEvt_Hole->Draw("box");
    h_IEtaIPhiILMapEvt_Hole->Draw("same,axis");
    h_IEtaIPhiILMapEvt_Ecut->SetMaximum(maxSeedEnergy);
    h_IEtaIPhiILMapEvt_Ecut->SetLineColor(kGray+1);
    h_IEtaIPhiILMapEvt_Ecut->Draw("box,same");
    
    drawLatexAdd(labelEnergy.Data(),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(labelDet,0.02,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    
  c3DBox->Print(Form("%s/3D/3D_Event%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));
  //************************************************************************************
  // draw 3D hist with box option   
  //************************************************************************************
    h_IEtaIPhiILMapEvt_Hole->Draw("box");
    h_IEtaIPhiILMapEvt_Hole->Draw("same,axis");
    h_IEtaIPhiILMapEvt_Ecut->Draw("box,same");
    for (int nCl = 0; nCl < nclusters; nCl++){
      if (h_IEtaIPhiILMapEvt_Cl[nCl]){
        h_IEtaIPhiILMapEvt_Cl[nCl]->SetMaximum(maxSeedEnergy);
        h_IEtaIPhiILMapEvt_Cl[nCl]->SetLineColor(colorCluster[nCl]);
        h_IEtaIPhiILMapEvt_Cl[nCl]->Draw("box,same");
        int cols = nCl%5;
        int rows = (nCl/5);
        drawLatexAdd(Form("%d: #it{E}_{cl} = %0.1f GeV", _mcpart_PDG[GetCorrectMCArrayEntry(rec_clusters.at(nCl).cluster_trueID)] , rec_clusters.at(nCl).cluster_E ), 0.02+0.2*cols, 0.97-(rows*0.55*textSizeLabelsRel), 0.55*textSizeLabelsRel, kFALSE, kFALSE, kFALSE, colorCluster[nCl]);
      }
    }
    drawLatexAdd(labelEnergy.Data(),0.98,0.02,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s, Event %d", labelDet.Data(), plotEvent),0.98,0.06,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (collisionsSys.Contains("Single")) 
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.02,0.015,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  
    c3DBox->Print(Form("%s/3DwClusters/3D_Event%d_Ecut_clusters.%s", outputDir.Data(),plotEvent, suffix.Data()));
  //************************************************************************************
  // draw 3D hist with box option   
  //************************************************************************************
    h_IEtaIPhiILMapEvt_Hole->Draw("box");
    h_IEtaIPhiILMapEvt_Hole->Draw("same,axis");
    h_IEtaIPhiILMapEvt_Ecut->Draw("box,same");
    for (int nMC = 0; nMC < nMCCurr; nMC++){
      if (h_IEtaIPhiILMapEvt_MCPart[nMC]){
        h_IEtaIPhiILMapEvt_MCPart[nMC]->SetMaximum(maxSeedEnergy);
        h_IEtaIPhiILMapEvt_MCPart[nMC]->SetLineColor(colorCluster[nMC]);
        h_IEtaIPhiILMapEvt_MCPart[nMC]->Draw("box,same");
        int cols = nMC%5;
        int rows = (nMC/5);
        drawLatexAdd((TString)(mcPartEnergies.at(nMC)).Data() , 0.02+0.2*cols, 0.98-(rows*0.42*textSizeLabelsRel), 0.42*textSizeLabelsRel, kFALSE, kFALSE, kFALSE, colorCluster[nMC]);
      }
    }
    drawLatexAdd(labelEnergy.Data(),0.98,0.02,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("%s, Event %d", labelDet.Data(), plotEvent),0.98,0.06,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (collisionsSys.Contains("Single")) 
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.02,0.015,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  
    c3DBox->Print(Form("%s/3DwClusters/3D_Event%d_Ecut_MCPart.%s", outputDir.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 1D projection on Z axis 
  //************************************************************************************
  TH1D* h_IZ_LFHCAL_ECut    = (TH1D*)h_IEtaIPhiILMapEvt_Ecut->Project3D("xe");
  TH1D* h_IZ_LFHCAL_Cl[50]  = {NULL};
  for (int nCl = 0; nCl < nclusters; nCl++){
    h_IZ_LFHCAL_Cl[nCl]     = (TH1D*)h_IEtaIPhiILMapEvt_Cl[nCl]->Project3D("xe");
  }
  DrawGammaCanvasSettings( c2DBox, 0.12, 0.01, 0.01, 0.1);
  c2DBox->cd();
  SetStyleHistoTH1ForGraphs(h_IZ_LFHCAL_ECut, "tower z index", "#it{E} (GeV)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.1, 505, 510);
    h_IZ_LFHCAL_ECut->Draw("hist,e1");
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IZ_LFHCAL_Cl[nCl]->SetLineColor(colorCluster[nCl]);
      h_IZ_LFHCAL_Cl[nCl]->Draw("hist,e1,same");
    }
    drawLatexAdd(Form("%s, %s",labelDet.Data(), labelEnergy.Data()),0.95,0.93,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
    if (collisionsSys.Contains("Single")){
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }
  
  c2DBox->Print(Form("%s/ZProjection_Event%d.%s", outputDir.Data(),plotEvent, suffix.Data()));


  TH2D* h_IEtaIPhi_LFHCAL_IZ[7]          = {NULL};
  TH2D* h_IEtaIPhi_LFHCAL_IZ_Cl[7][50]   = {{NULL}};
  
  for (Int_t z = 0; z < towersz; z++){
    h_IEtaIPhiILMapEvt_Ecut->GetXaxis()->SetRange(z+1,z+1);
    h_IEtaIPhi_LFHCAL_IZ[z]   = (TH2D*)h_IEtaIPhiILMapEvt_Ecut->Project3D("zye");
    h_IEtaIPhi_LFHCAL_IZ[z]->SetName(Form("etaphi_iL%d", z));
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IEtaIPhiILMapEvt_Cl[nCl]->GetXaxis()->SetRange(z+1,z+1);
      h_IEtaIPhi_LFHCAL_IZ_Cl[z][nCl]   = (TH2D*)h_IEtaIPhiILMapEvt_Cl[nCl]->Project3D(Form("cl%d_iL%d_zy_e", nCl, z));
    }
  }


  TCanvas* cPNG = new TCanvas("cPNG", "", (950*towersz+50)*4,950*4);
	split_canvas(cPNG, "cPNG1", towersz, true);
  Int_t padnum =0;

  for(Int_t z=0; z<towersz; z++){
    cPNG->cd(padnum+1);
    DrawVirtualPadSettings( cPNG->cd(padnum+1),  0.11, 0.03, 0.01, 0.1);


    h_IEtaIPhiMapEvt_Hole->Draw("box");
    h_IEtaIPhiMapEvt_Hole->Draw("axis,same");
    
    h_IEtaIPhi_LFHCAL_IZ[z]->SetLineColor(kGray+2);
    h_IEtaIPhi_LFHCAL_IZ[z]->Draw("box,same");

     for (int nCl = 0; nCl < nclusters; nCl++){
      h_IEtaIPhi_LFHCAL_IZ_Cl[z][nCl]->SetLineColor(colorCluster[nCl]);
      h_IEtaIPhi_LFHCAL_IZ_Cl[z][nCl]->Draw("box,same");
    }
   
    drawLatexAdd(Form("layer %d",z),0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
    if (z==towersz-1){
      drawLatexAdd(Form("%s, %s",labelDet.Data(),labelEnergy.Data()),0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
      if (collisionsSys.Contains("Single")){
        drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.14,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
    }

    padnum++;
  }
  cPNG->Print(Form("%s/3DwClusters/EtaPhi_Projection_Event%d.%s", outputDir.Data(), plotEvent, suffix.Data()));
  
  
}

