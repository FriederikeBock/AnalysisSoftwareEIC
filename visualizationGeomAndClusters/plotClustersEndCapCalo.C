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


void plotClustersEndCapCalo(
    TString nameFile          = "",
    TString inFileGeometry    = "",
    TString suffix            = "pdf",
    TString collisionsSys     = "PythiaMB",
    double energy             = -1.0,
    int plotEvent             = 7,
    int caloID                = 1,
    int clusterizerEnum       = 0,
    TString addName           = ""           
){
  
  
    //dimensions
  int towersx, towersy, offsetX, innerR;
  TString labelDet  = "";
  TString caloName = "";
  switch( caloID ){
    case 0:           // FHCAL
      towersx = 53; 
      towersy = 53;
      offsetX = 1;
      innerR  = 2;
      labelDet = Form("FHCAL%s, event %d", addName.Data(), plotEvent);
      caloName = "FHCAL";
      break;
    case 1:           // FEMC
      towersx = 198; 
      towersy = 198; 
      offsetX = 3;
      innerR  = 6;
      labelDet = Form("FEMC%s, event %d", addName.Data(), plotEvent);
      caloName = "FEMC";
      break;
    case 3:          // EEMC 
      towersx = 60; 
      towersy = 60; 
      offsetX = 0;
      innerR  = 4;
      labelDet = Form("EEMC%s, event %d", addName.Data(), plotEvent);
      caloName = "EEMC";
      break;
    case 5:          // EHCAL 
      towersx = 53; 
      towersy = 53; 
      offsetX = 0;
      innerR  = 2;
      labelDet = Form("EHCAL%s, event %d", addName.Data(), plotEvent);
      caloName = "EHCAL";
      break;
    case 9:          // EEMCG 
      towersx = 30; 
      towersy = 30; 
      offsetX = 0;
      innerR  = 2;
      labelDet = Form("EEMCG%s, event %d", addName.Data(), plotEvent);
      caloName = "EEMCG";
      break;
    default:
      std::cout << "Wrong caloID chosen" << std::endl;
      return;
  }

  int maxTower      = 0;
  if (towersx > towersy) 
    maxTower = towersx;
  else 
    maxTower = towersy;
  float minECutOff  = aggE[caloID]/2.;
  float seedEC      = seedE[caloID];
  float aggEC       = aggE[caloID];
  int verbosity     = 1;

  // dimensions 
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  
  Color_t colorCluster[50] = {kRed, kGreen+2, kBlue, kOrange, kCyan+1, kMagenta+1, kViolet-4, kTeal, kSpring-1, kYellow,
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
  } else if (collisionsSys.CompareTo("SingleE") == 0){
    labelEnergy   = "e^{-} simulation";
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
  TString outputDir                             = Form("plots/%s_%s_%s", dateForOutput.Data(), collisionsSys.Data(), caloName.Data());
  if (collisionsSys.Contains("Single"))
    outputDir = Form("%s_%.0f", outputDir.Data(), energy);
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDir+"/2D");
  gSystem->Exec("mkdir -p "+outputDir+"/2DwClusters");

  double textSizeSinglePad                    = 0.05;
  double textSizeLabelsPixel                  = 35;
  double textSizeLabelsRel                    = 58./1300;

  int granularity = 1;
  TH2F* h_IEtaIPhiMapEvt                   = NULL;
  TH2F* h_IEtaIPhiMapEvt_Ecut              = NULL;
  TH2F* h_IEtaIPhiMapEvt_Hole              = new TH2F("h_IEtaIPhiMapEvt_Hole", "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
  TH2F* h_IEtaIPhiMapEvt_Cl[50]            = {NULL};
  TH2F* h_IEtaIPhiMapEvt_MCPart[50]        = {NULL};
    
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
  geomStrct   arr_GeoStrt[towersx][towersy];  // array[iEta-1][iPhi-1]
  for(int i=0; i<(int)_calogeom_towers_N; i++){
    int tmpEta  = _calogeom_towers_iEta[i]-1;
    int tmpPhi  = _calogeom_towers_iPhi[i]-1;
    h_IEtaIPhiMapEvt_Hole->SetBinContent(tmpEta+2,tmpPhi+2, 1);
    arr_GeoStrt[tmpEta][tmpPhi].tower_Eta   = _calogeom_towers_Eta[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_Phi   = _calogeom_towers_Phi[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_x     = _calogeom_towers_x[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_y     = _calogeom_towers_y[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_z     = _calogeom_towers_z[i];
  }

  for (int x = 1; x < h_IEtaIPhiMapEvt_Hole->GetNbinsX()+1; x++){
    for (int y = 1; y < h_IEtaIPhiMapEvt_Hole->GetNbinsY()+1; y++){
      if (h_IEtaIPhiMapEvt_Hole->GetBinContent(x,y) > 0 )
        h_IEtaIPhiMapEvt_Hole->SetBinContent(x,y, 0);
      else 
        h_IEtaIPhiMapEvt_Hole->SetBinContent(x,y, 0.1);
    }
  }

  
  
  h_IEtaIPhiMapEvt       = new TH2F("h_IEtaIPhiMapEvt", "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
  h_IEtaIPhiMapEvt->Sumw2();
  h_IEtaIPhiMapEvt_Ecut  = new TH2F("h_IEtaIPhiMapEvt_Ecut", "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
  h_IEtaIPhiMapEvt_Ecut->Sumw2();
  for (int i = 0; i < 50; i++){
    h_IEtaIPhiMapEvt_Cl[i]      = new TH2F(Form("h_IEtaIPhiMapEvt_Cl%d",i), "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
    h_IEtaIPhiMapEvt_Cl[i]->Sumw2();
    h_IEtaIPhiMapEvt_MCPart[i]  = new TH2F(Form("h_IEtaIPhiMapEvt_MCPart%d",i), "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
    h_IEtaIPhiMapEvt_MCPart[i]->Sumw2();
  }
  
  
  float     sumedE          = 0;
  float     _nTowers        = 0;
  float*    _tower_E        = new float[_maxNTowersDR];
  int*      _tower_iEta       = new int[_maxNTowersDR];
  int*      _tower_iPhi       = new int[_maxNTowersDR];
  int*      _tower_trueID     = new int[_maxNTowersDR];

  Int_t nentries = Int_t(tt_event->GetEntries());
  for(Int_t i=0; i<nentries; i++){

    if( tt_event->LoadTree(i) < 0)  break;
    if( i!=plotEvent ) continue;
    tt_event->GetEntry(i);

    for(Int_t imc=0; (Int_t)imc<_nMCPart; imc++){
      TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
      _mcpart_Eta[imc]=truevec.Eta();
      _mcpart_Phi[imc]=truevec.Phi();
    }

    switch( caloID ){
      case 0:           // FHCAL
        _nTowers    = _nTowers_FHCAL;
        _tower_E    = _tower_FHCAL_E;
        _tower_iEta = _tower_FHCAL_iEta;
        _tower_iPhi = _tower_FHCAL_iPhi;
        _tower_trueID = _tower_FHCAL_trueID;
        break;
      case 1:           // FEMC
        _nTowers    = _nTowers_FEMC;
        _tower_E    = _tower_FEMC_E;
        _tower_iEta = _tower_FEMC_iEta;
        _tower_iPhi = _tower_FEMC_iPhi;
        _tower_trueID = _tower_FEMC_trueID;
        break;
      case 3:          // EEMC
        _nTowers    = _nTowers_EEMC;
        _tower_E    = _tower_EEMC_E;
        _tower_iEta = _tower_EEMC_iEta;
        _tower_iPhi = _tower_EEMC_iPhi;
        _tower_trueID = _tower_EEMC_trueID;
        break;
      case 5:          // EHCAL
        _nTowers    = _nTowers_EHCAL;
        _tower_E    = _tower_EHCAL_E;
        _tower_iEta = _tower_EHCAL_iEta;
        _tower_iPhi = _tower_EHCAL_iPhi;
        _tower_trueID = _tower_EHCAL_trueID;
        break;
      case 9:          // EEMCG
        _nTowers    = _nTowers_EEMCG;
        _tower_E    = _tower_EEMCG_E;
        _tower_iEta = _tower_EEMCG_iEta;
        _tower_iPhi = _tower_EEMCG_iPhi;
        _tower_trueID = _tower_EEMCG_trueID;
        break;
      default:
        std::cout << "Wrong caloID chosen" << std::endl;
        return;
      }

    if(_nTowers == 0 ) return;
  
    for(int itwr=0;(int)itwr<_nTowers;itwr++){
      sumedE += _tower_E[itwr];
      if (_tower_E[itwr] > 1e-3 && verbosity > 1) cout << _tower_iEta[itwr] << "\t" << _tower_iPhi[itwr] << "\t"<< _tower_E[itwr] << endl;
      h_IEtaIPhiMapEvt->Fill(_tower_iEta[itwr],_tower_iPhi[itwr],_tower_E[itwr]);
      if(_tower_E[itwr]>minECutOff){
        h_IEtaIPhiMapEvt_Ecut->Fill(_tower_iEta[itwr],_tower_iPhi[itwr],_tower_E[itwr]);
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
  
  // fill vector with towers for clusterization above aggregation threshold
  for(int itow=0; itow<(int)_nTowers; itow++){
    if(_tower_E[itow]>aggEC){
      towersStrct tempstructT;
      tempstructT.tower_E       = _tower_E[itow];
      tempstructT.tower_iEta    = _tower_iEta[itow];
      tempstructT.tower_iPhi    = _tower_iPhi[itow];
      tempstructT.tower_trueID  = _tower_trueID[itow];
      input_towers.push_back(tempstructT);
      input_towersForMC.push_back(tempstructT);
      if (verbosity > 1 ) std::cout << tempstructT.tower_E << "\t" << tempstructT.tower_iEta << "\t" << tempstructT.tower_iPhi << std::endl;
    }
  }

  // sort vector in descending energy order
  std::sort(input_towers.begin(), input_towers.end(), &acompare);
  std::sort(input_towersForMC.begin(), input_towersForMC.end(), &acompareTrueID);
  
  std::vector<int> clslabels;
  while (!input_towers.empty()) {
    cluster_towers.clear();
    clslabels.clear();

    clustersStrct tempstructC;

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

      cout << "filling hist for cluster: "  << nclusters << "\t E: " << tempstructC.cluster_E << "\t seed E: " << tempstructC.cluster_seed<<  endl;

      for (int tic = 0; tic < (int)cluster_towers.size(); tic++){
        int iEtaTwr = cluster_towers.at(tic).tower_iEta;
        int iPhiTwr = cluster_towers.at(tic).tower_iPhi;
        float eTwr  = cluster_towers.at(tic).tower_E;
        if (verbosity > 0) std::cout << "\t" << Form("%.3f",cluster_towers.at(tic).tower_E) << "\t eta \t" << cluster_towers.at(tic).tower_iEta << "\t phi \t" << cluster_towers.at(tic).tower_iPhi << std::endl;
        h_IEtaIPhiMapEvt_Cl[nclusters]->Fill(iEtaTwr, iPhiTwr, eTwr);
         
      }
      rec_clusters.push_back(tempstructC);
      recclusterEnergies.push_back(tempstructC.cluster_E);
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

  cout << "here" << endl;
  int nMCCurr   = 0;
  int nMCCurrID = input_towersForMC.at(0).tower_trueID;
  TString tempMCLabel = Form("%i: %.1f GeV, #eta = %1.2f", _mcpart_PDG[GetCorrectMCArrayEntry(nMCCurrID)] , _mcpart_E[GetCorrectMCArrayEntry(nMCCurrID)],_mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)]  );
  mcPartEnergies.push_back(tempMCLabel);
  while (!input_towersForMC.empty() && nMCCurr < 49) {
    if (nMCCurrID != input_towersForMC.at(0).tower_trueID){
      nMCCurr++;
      nMCCurrID = input_towersForMC.at(0).tower_trueID;
      tempMCLabel = Form("%i: %.1f GeV, #eta = %1.2f ", _mcpart_PDG[GetCorrectMCArrayEntry(nMCCurrID)] , _mcpart_E[GetCorrectMCArrayEntry(nMCCurrID)], _mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)] );
      mcPartEnergies.push_back(tempMCLabel);
    }

    int iEtaTwr = input_towersForMC.at(0).tower_iEta;
    int iPhiTwr = input_towersForMC.at(0).tower_iPhi;
    float eTwr    = input_towersForMC.at(0).tower_E;
    h_IEtaIPhiMapEvt_MCPart[nMCCurr]->Fill(iEtaTwr, iPhiTwr, eTwr);
    input_towersForMC.erase(input_towersForMC.begin());
  }
  
  
//   if (collisionsSys.Contains("Single"))
//     labelEnergy=Form("%s, #eta = %1.2f", labelEnergy.Data(), _mcpart_Eta[0]);
//   
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
      drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", recclusterEnergies.at(nCl) ), 0.14+0.2*cols, 0.93-(rows*0.55*textSizeLabelsRel), 0.55*textSizeLabelsRel, kFALSE, kFALSE, kFALSE, colorCluster[nCl]);
      
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
  
  for (int nMC = 0; nMC < nMCCurr+1; nMC++){
    if (h_IEtaIPhiMapEvt_MCPart[nMC]){
      h_IEtaIPhiMapEvt_MCPart[nMC]->SetLineColor(colorCluster[nMC]);
      h_IEtaIPhiMapEvt_MCPart[nMC]->Draw("box,same");
      int cols = nMC%4;
      int rows = (nMC/4);
      drawLatexAdd((TString)(mcPartEnergies.at(nMC)).Data() , 0.14+0.16*cols, 0.93-(rows*0.55*textSizeLabelsRel), 0.55*textSizeLabelsRel, kFALSE, kFALSE, kFALSE, colorCluster[nMC]);
    }
  }
  
  drawLatexAdd(labelDet.Data(),0.95,0.14+textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(labelEnergy.Data(),0.95,0.14,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  
  c2DBox->Print(Form("%s/2DwClusters/2D_Event%d_Ecut_MCPart.%s", outputDir.Data(),plotEvent, suffix.Data()));

  
}

