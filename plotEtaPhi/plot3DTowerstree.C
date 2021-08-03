#include <algorithm>
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
#include "TColor.h"

typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_iL;
} towersStrct;

typedef struct {
  float cluster_E;
  float cluster_seed;
  float cluster_Eta;
  float cluster_Phi;
  float cluster_Z;
  float cluster_X;
  float cluster_Y;
  float cluster_M02;
  float cluster_M20;
  int cluster_NTowers;
} clusterStrct;

// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }


void plot3DTowerstree(
    TString nameFile          = "",
    TString suffix            = "pdf",
    TString collisionsSys     = "PythiaMB",
    double energy             = -1.0,
    int plotEvent             = 7
){
  // dimensions 
  int towersx       = 524/5;
  int innerR        = 15/5;
  int offsetX       = 2;
  int towersz       = 5;
  float minECutOff  = 0.001;
  float seedE       = 0.25;
  float aggE        = 0.005;
  int verbosity     = 0;
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
  } else if (collisionsSys.CompareTo("SinglePion") == 0){
    labelEnergy   = Form("#pi^{-}, %0.1f GeV", energy);
  } else if (collisionsSys.CompareTo("SingleProton") == 0){
    labelEnergy   = Form("p, %0.1f GeV", energy);
  }
  
  TString dateForOutput                         = ReturnDateStringForOutput();

  TString outputDir                             = Form("plots/%s_%s",dateForOutput.Data(),collisionsSys.Data());
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
  TH2F* h_IEtaIPhiMapEvt_LFHCAL                   = NULL;
  TH2F* h_IEtaIPhiMapEvt_LFHCAL_Ecut              = NULL;
  TH2F* h_IEtaIPhiMapEvt_Hole                     = new TH2F("h_IEtaIPhiMapEvt_Hole", "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
  TH2F* h_IEtaIPhiMapEvt_LFHCAL_Cl[50]            = {NULL};
  TH3F* h_IEtaIPhiILMapEvt_LFHCAL                 = NULL;
  TH3F* h_IEtaIPhiILMapEvt_Hole                   = new TH3F("h_IEtaIPhiILMapEvt_Hole", "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5  );
  TH3F* h_IEtaIPhiILMapEvt_LFHCAL_Ecut            = NULL;
  TH3F* h_IEtaIPhiILMapEvt_LFHCAL_Cl[50]          = {NULL};
  
  for (int x = 0; x < towersx+1; x++){
    for (int y = 0; y < towersx+1; y++){
      if (TMath::Sqrt(TMath::Power(x-(towersx/2),2)+TMath::Power(y-(towersx/2),2)) > towersx/2){
        h_IEtaIPhiMapEvt_Hole->Fill(x,y, 0.1);
        for (int z = 0; z < 5; z++){
          h_IEtaIPhiILMapEvt_Hole->Fill(z,x,y, 0.1);
        }
      }
      if (TMath::Sqrt(TMath::Power(x-(double)((towersx+1+offsetX)/2),2)+TMath::Power(y-(double)((towersx+1)/2),2)) < (double)innerR){
        h_IEtaIPhiMapEvt_Hole->Fill(x,y, 0.1);
        for (int z = 0; z < 5; z++){
          h_IEtaIPhiILMapEvt_Hole->Fill(z,x,y, 0.1);
        }
      }
    }
  }
  
  // read file
  TChain* t_towers_LFHCAL = new TChain("event_tree", "event_tree");
  t_towers_LFHCAL->Add(nameFile.Data());
  
  const int _maxNTowersDR = towersx * towersx;
  int _nTowers_LFHCAL;
  float* _tower_LFHCAL_E           = new float[_maxNTowersDR];
  int* _tower_LFHCAL_iEta         = new int[_maxNTowersDR];
  int* _tower_LFHCAL_iPhi         = new int[_maxNTowersDR];
  int* _tower_LFHCAL_iL           = new int[_maxNTowersDR];
  int* _tower_LFHCAL_trueID       = new int[_maxNTowersDR];

  float sumedE                    = 0;
  if(t_towers_LFHCAL){
    h_IEtaIPhiMapEvt_LFHCAL       = new TH2F("h_IEtaIPhiMapEvt_LFHCAL", "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
    h_IEtaIPhiMapEvt_LFHCAL->Sumw2();
    h_IEtaIPhiMapEvt_LFHCAL_Ecut  = new TH2F("h_IEtaIPhiMapEvt_LFHCAL_Ecut", "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
    h_IEtaIPhiMapEvt_LFHCAL_Ecut->Sumw2();
    h_IEtaIPhiILMapEvt_LFHCAL       = new TH3F("h_IEtaIPhiILMapEvt_LFHCAL", "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5  );
    h_IEtaIPhiILMapEvt_LFHCAL->Sumw2();
    h_IEtaIPhiILMapEvt_LFHCAL_Ecut  = new TH3F("h_IEtaIPhiILMapEvt_LFHCAL_Ecut", "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
    h_IEtaIPhiILMapEvt_LFHCAL_Ecut->Sumw2();
    for (int i = 0; i < 50; i++){
      h_IEtaIPhiMapEvt_LFHCAL_Cl[i]  = new TH2F(Form("h_IEtaIPhiMapEvt_LFHCAL_Cl%d",i), "",  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
      h_IEtaIPhiMapEvt_LFHCAL_Cl[i]->Sumw2();
      h_IEtaIPhiILMapEvt_LFHCAL_Cl[i]  = new TH3F(Form("h_IEtaIPhiILMapEvt_LFHCAL_Cl%d",i), "", towersz, -0.5, towersz-0.5,  towersx+1, -0.5, towersx+0.5, towersx+1, -0.5, towersx+0.5);
      h_IEtaIPhiILMapEvt_LFHCAL_Cl[i]->Sumw2();
    }
    t_towers_LFHCAL->SetBranchAddress("tower_LFHCAL_N",             &_nTowers_LFHCAL);
    t_towers_LFHCAL->SetBranchAddress("tower_LFHCAL_E",             _tower_LFHCAL_E);
    t_towers_LFHCAL->SetBranchAddress("tower_LFHCAL_iEta",          _tower_LFHCAL_iEta);
    t_towers_LFHCAL->SetBranchAddress("tower_LFHCAL_iPhi",          _tower_LFHCAL_iPhi);
    t_towers_LFHCAL->SetBranchAddress("tower_LFHCAL_iL",            _tower_LFHCAL_iL);
    
    Int_t nentries = Int_t(t_towers_LFHCAL->GetEntries());
    for (Int_t i = 0; i < nentries; ++i) {
      if (t_towers_LFHCAL->LoadTree(i) < 0)
        break;
      if(i!=plotEvent) continue;
      t_towers_LFHCAL->GetEntry(i);
      for(int itwr=0;itwr<_nTowers_LFHCAL;itwr++){
        sumedE += _tower_LFHCAL_E[itwr];
        if (_tower_LFHCAL_E[itwr] > 1e-3 && verbosity) cout << _tower_LFHCAL_iEta[itwr] << "\t" << _tower_LFHCAL_iPhi[itwr] << "\t"<< _tower_LFHCAL_iL[itwr] << "\t"<< _tower_LFHCAL_E[itwr] << endl;
        h_IEtaIPhiMapEvt_LFHCAL->Fill(_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr],_tower_LFHCAL_E[itwr]);
        h_IEtaIPhiILMapEvt_LFHCAL->Fill(_tower_LFHCAL_iL[itwr],_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr], _tower_LFHCAL_E[itwr]);
        if(_tower_LFHCAL_E[itwr]>0.001){
          h_IEtaIPhiMapEvt_LFHCAL_Ecut->Fill(_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr],_tower_LFHCAL_E[itwr]);
          h_IEtaIPhiILMapEvt_LFHCAL_Ecut->Fill(_tower_LFHCAL_iL[itwr],_tower_LFHCAL_iEta[itwr],_tower_LFHCAL_iPhi[itwr], _tower_LFHCAL_E[itwr]);
        }  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
      }
    }
  }

  int nclusters = 0;

  // vector of towers used in the clusterizer
  std::vector<towersStrct> input_towers;
  // vector of towers within the currently found cluster
  std::vector<towersStrct> cluster_towers;
  std::vector<clusterStrct> rec_clusters;
  
  // fill vector with towers for clusterization above aggregation threshold
  for(int itow=0; itow<_nTowers_LFHCAL; itow++){
    if(_tower_LFHCAL_E[itow]>aggE){
      towersStrct tempstructT;
      tempstructT.tower_E       = _tower_LFHCAL_E[itow];
      tempstructT.tower_iEta    = _tower_LFHCAL_iEta[itow];
      tempstructT.tower_iPhi    = _tower_LFHCAL_iPhi[itow];
      tempstructT.tower_iL      = _tower_LFHCAL_iL[itow];
      input_towers.push_back(tempstructT);
      if (verbosity) std::cout << tempstructT.tower_E << "\t" << tempstructT.tower_iEta << "\t" << tempstructT.tower_iPhi << "\t" << tempstructT.tower_iL << std::endl;
    }
  }

  // sort vector in descending energy order
  std::sort(input_towers.begin(), input_towers.end(), &acompare);
  while (!input_towers.empty()) {
    cluster_towers.clear();
    clusterStrct tempstructC;
    tempstructC.cluster_E       = 0;
    tempstructC.cluster_seed    = 0;
    tempstructC.cluster_Eta     = -1000;
    tempstructC.cluster_Phi     = -1000;
    tempstructC.cluster_Z       = -1000;
    tempstructC.cluster_X       = -1000;
    tempstructC.cluster_Y       = -1000;
    tempstructC.cluster_M02     = -1000;
    tempstructC.cluster_M20     = -1000;
    tempstructC.cluster_NTowers = 0;

    // always start with highest energetic tower
    if(input_towers.at(0).tower_E > seedE){
      if (verbosity) std::cout << "new cluster" << std::endl;
      // fill seed cell information into current cluster
      tempstructC.cluster_E       = input_towers.at(0).tower_E;
      tempstructC.cluster_seed    = tempstructC.cluster_E;
      tempstructC.cluster_NTowers = 1;
      cluster_towers.push_back(input_towers.at(0));
      if (verbosity) std::cout << "seed: "<<  input_towers.at(0).tower_iEta << "\t" << input_towers.at(0).tower_iPhi << "\t" << input_towers.at(0).tower_iL << "\t E:"<< tempstructC.cluster_E << std::endl;
      

      // remove seed tower from sample
      input_towers.erase(input_towers.begin());
      for (int tit = 0; tit < cluster_towers.size(); tit++){
        // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
        int refC = 0;
        int iEtaTwr = cluster_towers.at(tit).tower_iEta;
        int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
        int iLTwr   = cluster_towers.at(tit).tower_iL;
        for (int ait = 0; ait < input_towers.size(); ait++){
          int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
          int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
          int iLTwrAgg   = input_towers.at(ait).tower_iL;
          // first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
          int deltaL    = TMath::Abs(iLTwrAgg-iLTwr) ;
          int deltaPhi  = TMath::Abs(iPhiTwrAgg-iPhiTwr) ;
          int deltaEta  = TMath::Abs(iEtaTwrAgg-iEtaTwr) ;
          bool neighbor = (deltaL+deltaPhi+deltaEta == 1);
          bool corner2D = (deltaL == 0 && deltaPhi == 1 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 0 && deltaEta == 1) || (deltaL == 1 && deltaPhi == 1 && deltaEta == 0);          
          if(neighbor || corner2D ){
            // only aggregate towers with lower energy than current tower
//             if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + (aggE/10.))) continue;
            tempstructC.cluster_E+=input_towers.at(ait).tower_E;
            tempstructC.cluster_NTowers++;
            cluster_towers.push_back(input_towers.at(ait));
            if (verbosity) std::cout << "aggregated: "<< iEtaTwrAgg << "\t" << iPhiTwrAgg << "\t" << iLTwrAgg << "\t E:" << input_towers.at(ait).tower_E << "\t reference: "<< refC << "\t"<< iEtaTwr << "\t" << iPhiTwr << "\t" << iLTwr << "\t cond.: \t"<< neighbor << "\t" << corner2D << "\t  diffs: " << deltaEta << "\t" << deltaPhi << "\t" << deltaL<< std::endl;
            input_towers.erase(input_towers.begin()+ait);
            ait--;
            refC++;
          }
        }
      }
      cout << "filling hist for cluster: "  << nclusters << "\t E: " << tempstructC.cluster_E << "\t seed E: " << tempstructC.cluster_seed<<  endl;
      for (int tic = 0; tic < cluster_towers.size(); tic++){
        int iEtaTwr = cluster_towers.at(tic).tower_iEta;
        int iPhiTwr = cluster_towers.at(tic).tower_iPhi;
        int iLTwr   = cluster_towers.at(tic).tower_iL;
        float eTwr    = cluster_towers.at(tic).tower_E;
        if (verbosity) cout << iEtaTwr << "\t"<< iPhiTwr << "\t"<< iLTwr << "\t" << eTwr << endl;
        h_IEtaIPhiMapEvt_LFHCAL_Cl[nclusters]->Fill(iEtaTwr, iPhiTwr, eTwr);
        h_IEtaIPhiILMapEvt_LFHCAL_Cl[nclusters]->Fill(iLTwr, iEtaTwr, iPhiTwr, eTwr);
      }
      rec_clusters.push_back(tempstructC);

      nclusters++;
      
    } else {
      input_towers.clear();
    }
  }
  
  float firstClusterE     = 0;
  if(nclusters==0) 
    firstClusterE   = 0;
  else 
    firstClusterE   = rec_clusters.at(0).cluster_E;
    
  float maxSeedEnergy = 0;
  for (int it = 0; it < rec_clusters.size(); it++){
      if (rec_clusters.at(it).cluster_seed > maxSeedEnergy)
        maxSeedEnergy   =  rec_clusters.at(it).cluster_seed;
  }
  
  if (sumedE < 0.01){
    cout << "total E calo: " << sumedE <<  " No particle hit the calorimeter! Aborting." << endl; 
    return;
  } else {
    cout << "total E calo: " << sumedE << endl; 
  }
  
  // 2D PLOT
  TCanvas* cReso = new TCanvas("cReso","",0,0,1000,950);
  DrawGammaCanvasSettings( cReso, 0.11, 0.1, 0.01, 0.1);

  // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_LFHCAL_Ecut, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
    
    h_IEtaIPhiMapEvt_LFHCAL_Ecut->Draw("colz");
    h_IEtaIPhiMapEvt_Hole->SetLineColor(kGray);
    h_IEtaIPhiMapEvt_Hole->Draw("box,same");
    
    drawLatexAdd(Form("LHCal, %s",labelEnergy.Data()),0.88,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/2D/2D_Event%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));

  // reset 2D canvas
  DrawGammaCanvasSettings( cReso, 0.11, 0.03, 0.01, 0.1);
  
  //************************************************************************************
  // draw 2D hist with box option
  //************************************************************************************  
    SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_Hole, "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
    h_IEtaIPhiMapEvt_Hole->GetZaxis()->SetRangeUser(0,h_IEtaIPhiMapEvt_LFHCAL_Ecut->GetMaximum()*1.1);
    h_IEtaIPhiMapEvt_Hole->Draw("box");
    h_IEtaIPhiMapEvt_Hole->Draw("same,axis");
    h_IEtaIPhiMapEvt_LFHCAL_Ecut->SetLineColor(kGray+2);
    h_IEtaIPhiMapEvt_LFHCAL_Ecut->Draw("box,same");
    
    
    drawLatexAdd(Form("LHCal, %s",labelEnergy.Data()),0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  

  cReso->Print(Form("%s/2D/2D_Event%d_Ecut_Box.%s", outputDir.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 2D hist with box option and different clusters on top
  //************************************************************************************
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IEtaIPhiMapEvt_LFHCAL_Cl[nCl]->SetLineColor(colorCluster[nCl]);
      h_IEtaIPhiMapEvt_LFHCAL_Cl[nCl]->Draw("box,same");
    }
    if (collisionsSys.Contains("Single")){
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.14,0.93,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      if(nclusters>0)drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", firstClusterE ),0.14,0.88,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
  cReso->Print(Form("%s/2DwClusters/2D_Event%d_Ecut_clusters.%s", outputDir.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 3D hist with box option   
  //************************************************************************************
  DrawGammaCanvasSettings( cReso, 0.12, 0.05, 0.01, 0.1);
  cReso->SetTheta(16);
  cReso->SetPhi(68);
   
    SetStyleHistoTH3ForGraphs(h_IEtaIPhiILMapEvt_Hole, "tower z index", "tower #phi index","tower #eta index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1, 1.1, 1.15, 505, 510,510);

    h_IEtaIPhiILMapEvt_Hole->SetLineColor(kGray);
    h_IEtaIPhiILMapEvt_Hole->Draw("box");
    h_IEtaIPhiILMapEvt_Hole->Draw("same,axis");
    h_IEtaIPhiILMapEvt_LFHCAL_Ecut->SetMaximum(maxSeedEnergy);
    h_IEtaIPhiILMapEvt_LFHCAL_Ecut->SetLineColor(kGray+1);
    h_IEtaIPhiILMapEvt_LFHCAL_Ecut->Draw("box,same");
    
    drawLatexAdd(labelEnergy.Data(),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd("LFHCal",0.02,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    
  cReso->Print(Form("%s/3D/3D_Event%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));

  //************************************************************************************
  // draw 3D hist with box option   
  //************************************************************************************
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IEtaIPhiILMapEvt_LFHCAL_Cl[nCl]->SetMaximum(maxSeedEnergy);
      h_IEtaIPhiILMapEvt_LFHCAL_Cl[nCl]->SetLineColor(colorCluster[nCl]);
      h_IEtaIPhiILMapEvt_LFHCAL_Cl[nCl]->Draw("box,same");
    }
    drawLatexAdd(labelEnergy.Data(),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd("LFHCal",0.02,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    if (collisionsSys.Contains("Single")) 
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV, #it{E}_{cl} = %0.1f GeV",sumedE,firstClusterE  ),0.02,0.015,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  
    cReso->Print(Form("%s/3DwClusters/3D_Event%d_Ecut_clusters.%s", outputDir.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 1D projection on Z axis 
  //************************************************************************************
  TH1D* h_IZ_LFHCAL_ECut    = (TH1D*)h_IEtaIPhiILMapEvt_LFHCAL_Ecut->Project3D("xe");
  TH1D* h_IZ_LFHCAL_Cl[50]  = {NULL};
  for (int nCl = 0; nCl < nclusters; nCl++){
    h_IZ_LFHCAL_Cl[nCl]     = (TH1D*)h_IEtaIPhiILMapEvt_LFHCAL_Cl[nCl]->Project3D("xe");
  }
  DrawGammaCanvasSettings( cReso, 0.12, 0.01, 0.01, 0.1);
  SetStyleHistoTH1ForGraphs(h_IZ_LFHCAL_ECut, "tower z index", "#it{E} (GeV)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9, 1.1, 505, 510);
    h_IZ_LFHCAL_ECut->Draw("hist,e1");
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IZ_LFHCAL_Cl[nCl]->SetLineColor(colorCluster[nCl]);
      h_IZ_LFHCAL_Cl[nCl]->Draw("hist,e1,same");
    }
    drawLatexAdd(Form("LHCal, %s",labelEnergy.Data()),0.95,0.93,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
    if (collisionsSys.Contains("Single")){
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.95,0.88,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      if(nclusters>0) drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", firstClusterE ),0.95,0.83,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }
  
  cReso->Print(Form("%s/ZProjection_Event%d.%s", outputDir.Data(),plotEvent, suffix.Data()));


  TH2D* h_IEtaIPhi_LFHCAL_IZ[5]          = {NULL};
  TH2D* h_IEtaIPhi_LFHCAL_IZ_Cl[5][50]   = {{NULL}};
  
  for (Int_t z = 0; z < 5; z++){
    h_IEtaIPhiILMapEvt_LFHCAL_Ecut->GetXaxis()->SetRange(z+1,z+1);
    h_IEtaIPhi_LFHCAL_IZ[z]   = (TH2D*)h_IEtaIPhiILMapEvt_LFHCAL_Ecut->Project3D("zye");
    h_IEtaIPhi_LFHCAL_IZ[z]->SetName(Form("etaphi_iL%d", z));
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IEtaIPhiILMapEvt_LFHCAL_Cl[nCl]->GetXaxis()->SetRange(z+1,z+1);
      h_IEtaIPhi_LFHCAL_IZ_Cl[z][nCl]   = (TH2D*)h_IEtaIPhiILMapEvt_LFHCAL_Cl[nCl]->Project3D(Form("cl%d_iL%d_zy_e", nCl, z));
    }
  }

  
  TCanvas* cPNG = new TCanvas("cPNG", "", (950*5+50)*4,950*4);
	split_canvas(cPNG, "cPNG1", 5, true);
  Int_t padnum =0;

  for(Int_t z=0; z<5; z++){
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
    if (z==4){
      drawLatexAdd(Form("LHCal, %s",labelEnergy.Data()),0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
      if (collisionsSys.Contains("Single")){
        drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.14,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        if(nclusters>0)drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", firstClusterE ),0.14,0.87,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      }
    }

    padnum++;
  }
  cPNG->Print(Form("%s/3DwClusters/EtaPhi_Projection_Event%d.%s", outputDir.Data(), plotEvent, suffix.Data()));
  
  
}

