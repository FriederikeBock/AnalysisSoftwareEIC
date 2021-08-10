#include <algorithm>
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
#include "TColor.h"
#include <TChain.h>
#include <TVector3.h>

typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
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


void plot2DHCalInClusters(
    TString nameFile          = "",
    TString suffix            = "pdf",
    TString collisionsSys     = "PythiaMB",
    double energy             = -1.0,
    int plotEvent             = 7,
    TString addName           = ""           
){
  // dimensions 
  int towersx       = 26;
  int towersy       = 64;
  int maxTower      = 0;
  if (towersx > towersy) 
    maxTower = towersx;
  else 
    maxTower = towersy;
  float minECutOff  = 0.001;
  float seedE       = 0.2;
  float aggE        = 0.05;
  int verbosity     = 1;
  TString labelDet  = Form("HCALin%s, event %d", addName.Data(), plotEvent);
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
  } else if (collisionsSys.CompareTo("SingleGamma") == 0){
    labelEnergy   = Form("#gamma, %0.1f GeV", energy);
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
  TString labelParticlProp = "";
  
  TString dateForOutput                         = ReturnDateStringForOutput();

  TString outputDir                             = Form("plots/HCALin%s_%s",dateForOutput.Data(),collisionsSys.Data());
  if (collisionsSys.Contains("Single"))
    outputDir = Form("%s_%.0f", outputDir.Data(), energy);
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDir+"/2D");
  gSystem->Exec("mkdir -p "+outputDir+"/2DwClusters");

  double textSizeSinglePad                    = 0.05;
  double textSizeLabelsPixel                  = 35;
  double textSizeLabelsRel                    = 58./1300;

  int granularity = 1;
  TH2F* h_IEtaIPhiMapEvt_HCALIN                   = NULL;
  TH2F* h_IEtaIPhiMapEvt_HCALIN_Ecut              = NULL;
  TH2F* h_IEtaIPhiMapEvt_HCALIN_Cl[50]            = {NULL};
    
  // read file
  TChain* t_towers_HCALIN = new TChain("event_tree", "event_tree");
  t_towers_HCALIN->Add(nameFile.Data());
  
  const int _maxNTowersDR = towersx * towersy;
  int _nTowers_HCALIN;
  float* _tower_HCALIN_E           = new float[_maxNTowersDR];
  int* _tower_HCALIN_iEta         = new int[_maxNTowersDR];
  int* _tower_HCALIN_iPhi         = new int[_maxNTowersDR];
  int* _tower_HCALIN_trueID       = new int[_maxNTowersDR];
  const int _maxNMCPart = 10000;
  int _nMCPart;
  float* _mcpart_E                = new float[_maxNMCPart];
  float* _mcpart_px                = new float[_maxNMCPart];
  float* _mcpart_py                = new float[_maxNMCPart];
  float* _mcpart_pz                = new float[_maxNMCPart];
  float* _mcpart_Eta              = new float[_maxNMCPart];
  float* _mcpart_Phi              = new float[_maxNMCPart];
  
  float sumedE                    = 0;
  if(t_towers_HCALIN){
    h_IEtaIPhiMapEvt_HCALIN       = new TH2F("h_IEtaIPhiMapEvt_HCALIN", "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
    h_IEtaIPhiMapEvt_HCALIN->Sumw2();
    h_IEtaIPhiMapEvt_HCALIN_Ecut  = new TH2F("h_IEtaIPhiMapEvt_HCALIN_Ecut", "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
    h_IEtaIPhiMapEvt_HCALIN_Ecut->Sumw2();
    for (int i = 0; i < 50; i++){
      h_IEtaIPhiMapEvt_HCALIN_Cl[i]  = new TH2F(Form("h_IEtaIPhiMapEvt_HCALIN_Cl%d",i), "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
      h_IEtaIPhiMapEvt_HCALIN_Cl[i]->Sumw2();
    }
    t_towers_HCALIN->SetBranchAddress("tower_HCALIN_N",             &_nTowers_HCALIN);
    t_towers_HCALIN->SetBranchAddress("tower_HCALIN_E",             _tower_HCALIN_E);
    t_towers_HCALIN->SetBranchAddress("tower_HCALIN_iEta",          _tower_HCALIN_iEta);
    t_towers_HCALIN->SetBranchAddress("tower_HCALIN_iPhi",          _tower_HCALIN_iPhi);
    t_towers_HCALIN->SetBranchAddress("nMCPart",       &_nMCPart);
    t_towers_HCALIN->SetBranchAddress("mcpart_px",     _mcpart_px);
    t_towers_HCALIN->SetBranchAddress("mcpart_py",     _mcpart_py);
    t_towers_HCALIN->SetBranchAddress("mcpart_pz",     _mcpart_pz);
    t_towers_HCALIN->SetBranchAddress("mcpart_E",      _mcpart_E);
    
    Int_t nentries = Int_t(t_towers_HCALIN->GetEntries());
    for (Int_t i = 0; i < nentries; ++i) {
      
      if (t_towers_HCALIN->LoadTree(i) < 0)
        break;
      if(i!=plotEvent) continue;
      t_towers_HCALIN->GetEntry(i);
      
      for(Int_t imc=0; (Int_t)imc<_nMCPart; imc++){
        TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
        _mcpart_Eta[imc]=truevec.Eta();
        _mcpart_Phi[imc]=truevec.Phi();
      }
      
      for(int itwr=0;(int)itwr<_nTowers_HCALIN;itwr++){
        sumedE += _tower_HCALIN_E[itwr];
        if (_tower_HCALIN_E[itwr] > 1e-3 && verbosity) cout << _tower_HCALIN_iEta[itwr] << "\t" << _tower_HCALIN_iPhi[itwr] << "\t"<< _tower_HCALIN_E[itwr] << endl;
        h_IEtaIPhiMapEvt_HCALIN->Fill(_tower_HCALIN_iEta[itwr],_tower_HCALIN_iPhi[itwr],_tower_HCALIN_E[itwr]);
        if(_tower_HCALIN_E[itwr]>0.001){
          h_IEtaIPhiMapEvt_HCALIN_Ecut->Fill(_tower_HCALIN_iEta[itwr],_tower_HCALIN_iPhi[itwr],_tower_HCALIN_E[itwr]);
        }
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
  for(int itow=0; itow<(int)_nTowers_HCALIN; itow++){
    if(_tower_HCALIN_E[itow]>aggE){
      towersStrct tempstructT;
      tempstructT.tower_E       = _tower_HCALIN_E[itow];
      tempstructT.tower_iEta    = _tower_HCALIN_iEta[itow];
      tempstructT.tower_iPhi    = _tower_HCALIN_iPhi[itow];
      input_towers.push_back(tempstructT);
      if (verbosity) std::cout << tempstructT.tower_E << "\t" << tempstructT.tower_iEta << "\t" << tempstructT.tower_iPhi << std::endl;
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
      if (verbosity) std::cout << "seed: "<<  input_towers.at(0).tower_iEta << "\t" << input_towers.at(0).tower_iPhi << "\t E:"<< tempstructC.cluster_E << std::endl;
      

      // remove seed tower from sample
      input_towers.erase(input_towers.begin());
      for (int tit = 0; tit < (int)cluster_towers.size(); tit++){
        // Now go recursively to all neighbours and add them to the cluster if they fulfill the conditions
        int refC = 0;
        int iEtaTwr = cluster_towers.at(tit).tower_iEta;
        int iPhiTwr = cluster_towers.at(tit).tower_iPhi;
        for (int ait = 0; ait < (int)input_towers.size(); ait++){
          int iEtaTwrAgg = input_towers.at(ait).tower_iEta;
          int iPhiTwrAgg = input_towers.at(ait).tower_iPhi;
          
          if (iPhiTwr < 5 && iPhiTwrAgg > towersy-5){
            iPhiTwrAgg= iPhiTwrAgg-towersy;
          }
          if (iPhiTwr > towersy-5 && iPhiTwrAgg < 5){
            iPhiTwr= iPhiTwr-towersy;
          }

          // first condition asks for V3-like neighbors, while second condition also checks diagonally attached towers
          int deltaPhi  = TMath::Abs(iPhiTwrAgg-iPhiTwr) ;
          int deltaEta  = TMath::Abs(iEtaTwrAgg-iEtaTwr) ;
          bool neighbor = (deltaPhi+deltaEta == 1);
          bool corner2D = deltaPhi == 1 && deltaEta == 1;
          if(neighbor || corner2D ){
            // only aggregate towers with lower energy than current tower
//             if(input_towers.at(ait).tower_E >= (cluster_towers.at(tit).tower_E + (aggE/10.))) continue;
            tempstructC.cluster_E+=input_towers.at(ait).tower_E;
            tempstructC.cluster_NTowers++;
            cluster_towers.push_back(input_towers.at(ait));
            if (verbosity) std::cout  << "aggregated: "<< iEtaTwrAgg << "\t" << iPhiTwrAgg << "\t E:" << input_towers.at(ait).tower_E 
                                      << "\t reference: "<< refC << "\t"<< iEtaTwr << "\t" << iPhiTwr 
                                      << "\t cond.: \t"<< neighbor << "\t" << corner2D << "\t  diffs: " << deltaEta << "\t" << deltaPhi << std::endl;
            input_towers.erase(input_towers.begin()+ait);
            ait--;
            refC++;
          }
        }
      }
      cout << "filling hist for cluster: "  << nclusters << "\t E: " << tempstructC.cluster_E << "\t seed E: " << tempstructC.cluster_seed<<  endl;
      for (int tic = 0; tic < (int)cluster_towers.size(); tic++){
        int iEtaTwr = cluster_towers.at(tic).tower_iEta;
        int iPhiTwr = cluster_towers.at(tic).tower_iPhi;
        float eTwr    = cluster_towers.at(tic).tower_E;
        if (verbosity) cout << iEtaTwr << "\t"<< iPhiTwr << "\t" << eTwr << endl;
        h_IEtaIPhiMapEvt_HCALIN_Cl[nclusters]->Fill(iEtaTwr, iPhiTwr, eTwr);
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
  for (int it = 0; it < (int)rec_clusters.size(); it++){
      if (rec_clusters.at(it).cluster_seed > maxSeedEnergy)
        maxSeedEnergy   =  rec_clusters.at(it).cluster_seed;
  }
  
  if (sumedE < 0.01){
    cout << "total E calo: " << sumedE <<  " No particle hit the calorimeter! Aborting." << endl; 
    return;
  } else {
    cout << "total E calo: " << sumedE << endl; 
  }
  
  if (collisionsSys.Contains("Single"))
    labelParticlProp=Form("#eta = %1.2f, #varphi = %1.3f", _mcpart_Eta[0], _mcpart_Phi[0]);
  
  // 2D PLOT

  // 2D PLOT
  Int_t pixelBaseCanvas = 1000;
  Float_t fracX = (Float_t)towersx/maxTower;
  Float_t fracY = (Float_t)towersy/maxTower;
  Int_t pixelX  = 35*3.4+fracX*pixelBaseCanvas;
  Int_t pixelY  = 35*1.6+fracY*pixelBaseCanvas;
  Float_t marginBottom = (Float_t)(35*1.5)/(Float_t)pixelY;
  Float_t marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
  Float_t marginLeft = 35*1.5/(Float_t)pixelX;
  Float_t marginRight = 35*1.9/(Float_t)pixelX;
  
  cout << pixelX << "\t" << pixelY << "\t" << marginLeft << "\t" << marginRight << "\t" << marginTop << "\t" << marginBottom << endl;
  
  TCanvas* c2DColz = new TCanvas("c2DColz","",0,0,pixelX,pixelY);
  DrawGammaCanvasSettings( c2DColz, marginLeft, marginRight, marginTop, marginBottom);
  c2DColz->SetLogz(1);

  // 1D PLOT
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_HCALIN_Ecut, "tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.45,1.05);
  h_IEtaIPhiMapEvt_HCALIN_Ecut->GetXaxis()->SetLabelOffset(-0.015);
  h_IEtaIPhiMapEvt_HCALIN_Ecut->GetXaxis()->SetTickLength(0.02);
  
    h_IEtaIPhiMapEvt_HCALIN_Ecut->Draw("colz");
    
    drawLatexAdd(Form("%s, %s",labelDet.Data(), labelEnergy.Data()),1-marginRight-0.05, 1-marginTop-0.035,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(labelParticlProp,1-marginRight-0.05,1-marginTop-0.035-0.6*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    if (collisionsSys.Contains("Single")){
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),1-marginRight-0.05,1-marginTop-0.035-2*0.6*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }

  c2DColz->Print(Form("%s/2D/2D_Event%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()));

  // reset 2D canvas
  pixelX  = 35*2.7+fracX*pixelBaseCanvas;
  pixelY  = 35*1.5+fracY*pixelBaseCanvas;
  marginBottom = (Float_t)(35*1.4)/(Float_t)pixelY;
  marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
  marginLeft = 35*1.4/(Float_t)pixelX;
  marginRight = 35*0.3/(Float_t)pixelX;
  TCanvas* c2DBox = new TCanvas("c2DBox","",0,0,pixelX,pixelY);

  DrawGammaCanvasSettings( c2DBox, marginLeft, marginRight, marginTop, marginBottom);
  
  //************************************************************************************
  // draw 2D hist with box option
  //************************************************************************************  
  SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_HCALIN_Ecut, "tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.4,1.0);
  
    h_IEtaIPhiMapEvt_HCALIN_Ecut->SetLineColor(kGray+2);
    h_IEtaIPhiMapEvt_HCALIN_Ecut->Draw("box");
    
    drawLatexAdd(Form("%s, %s",labelDet.Data(), labelEnergy.Data()),1-marginRight-0.05,marginBottom+0.03+0.5*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
    drawLatexAdd(labelParticlProp,1-marginRight-0.05,marginBottom+0.03,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
    if (collisionsSys.Contains("Single")){
      drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),marginLeft+0.04,1-marginTop-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
  c2DBox->Print(Form("%s/2D/2D_Event%d_Ecut_Box.%s", outputDir.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 2D hist with box option and different clusters on top
  //************************************************************************************
    for (int nCl = 0; nCl < nclusters; nCl++){
      h_IEtaIPhiMapEvt_HCALIN_Cl[nCl]->SetLineColor(colorCluster[nCl]);
      h_IEtaIPhiMapEvt_HCALIN_Cl[nCl]->Draw("box,same");
    }
    if (collisionsSys.Contains("Single")){
      if(nclusters>0)drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", firstClusterE ),marginLeft+0.04,1-marginTop-2*0.85*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
  c2DBox->Print(Form("%s/2DwClusters/2D_Event%d_Ecut_clusters.%s", outputDir.Data(),plotEvent, suffix.Data()));
    
}

