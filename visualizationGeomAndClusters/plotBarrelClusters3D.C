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

float GetZ(float eta, float X, float Y, float Z){ // X,Y,Z - center bin
  float theta = 2 * TMath::ATan( TMath::Exp(-eta) );
  float R = TMath::Sqrt( X*X + Y*Y + Z*Z);
  return R * TMath::Cos(theta);
}

void plotBarrelClusters3D(
    TString nameFile          = "/home/ewa/alice/EIC/G4EICDetector_eventtree.root",
    TString inFileGeometry    = "/home/ewa/alice/EIC/geometry_ideal_BECALproj_EEMC_LFHCAL_corrected.root",
    TString suffix            = "pdf",
    TString collisionsSys     = "PythiaMB",
    double energy             = -1.0,
    int plotEvent             = 0,
    int caloID                = 6,     // 6 - HCAL in, 7 - HCAL out, 10 - BECAL
    int clusterizerEnum       = 0,
    TString addName           = ""           
){
  //dimensions
  int towersx, towersy;
  TString labelDet  = "";
  TString caloName = "";
  switch( caloID ){
    case 6:           // HCAL in
      towersx = 26; 
      towersy = 64;
      labelDet = Form("HCALin%s, event %d", addName.Data(), plotEvent);
      caloName = "HCALin";
      break;
    case 7:           // HCAL out
      towersx = 24; 
      towersy = 64; 
      labelDet = Form("HCALout%s, event %d", addName.Data(), plotEvent);
      caloName = "HCALout";
      break;
    case 10:          // BECAL 
      towersx = 70; 
      towersy = 128; 
      labelDet = Form("BECAL%s, event %d", addName.Data(), plotEvent);
      caloName = "BECAL";
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
  float seedEC       = seedE[caloID];
  float aggEC        = aggE[caloID];
  int verbosity     = 1;

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
  TString outputDirLocal                        = Form("plots/%s_%s_%s",dateForOutput.Data(),collisionsSys.Data(),caloName.Data());
  std::cout << outputDirLocal.Data() << std::endl;
  if (collisionsSys.Contains("Single"))
    outputDirLocal = Form("%s_%.0f", outputDirLocal.Data(), energy);
  gSystem->Exec("mkdir -p "+outputDirLocal);
  gSystem->Exec("mkdir -p "+outputDirLocal+"/2D");
  gSystem->Exec("mkdir -p "+outputDirLocal+"/2DwClusters");
  gSystem->Exec("mkdir -p "+outputDirLocal+"/3D");

  double textSizeSinglePad                    = 0.05;
  double textSizeLabelsPixel                  = 35;
  double textSizeLabelsRel                    = 58./1300;

  int granularity = 1;
  TH2F* h_IEtaIPhiMapEvt                      = NULL;
  TH2F* h_IEtaIPhiMapEvt_Ecut                 = NULL;
  TH2F* h_IEtaIPhiMapEvt_Cl[50]               = {NULL};
  TH2F* h_IEtaIPhiMapEvt_MCPart[50]           = {NULL};

  float* _mcpart_Phi                          = new float[_maxNMCPart];

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
    arr_GeoStrt[tmpEta][tmpPhi].tower_Eta   = _calogeom_towers_Eta[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_Phi   = _calogeom_towers_Phi[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_x     = _calogeom_towers_x[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_y     = _calogeom_towers_y[i];
    arr_GeoStrt[tmpEta][tmpPhi].tower_z     = _calogeom_towers_z[i];
  }

  // create the binning for the z component
  std::vector<Double_t> vec_zBins(towersx+1,0);
  float DeltaEta  = 1.0 * TMath::Abs(_calogeom_towers_Eta[0] - _calogeom_towers_Eta[_calogeom_towers_N-1] ) / towersx;
  float EtaMax;

  if(caloID==10)    EtaMax    = _calogeom_towers_Eta[0] + 0.5*DeltaEta;
  else  EtaMax = _calogeom_towers_Eta[_calogeom_towers_N-1] + 0.5*DeltaEta;

  for(int i=0; i<towersx; i++){
    vec_zBins.at(i) = GetZ( EtaMax - i*DeltaEta, _calogeom_towers_x[i*towersy], _calogeom_towers_y[i*towersy], _calogeom_towers_z[i*towersy]);
  }
  vec_zBins.at(towersx) =  vec_zBins[towersx-1] - TMath::Abs(vec_zBins[towersx-1] - vec_zBins[towersx-2]);   // approx for the last bin
  std::sort( vec_zBins.begin(), vec_zBins.end() );

  Double_t zBins[vec_zBins.size()];
  std::copy(vec_zBins.begin(), vec_zBins.end(), zBins);

  //defining the histograms
  h_IEtaIPhiMapEvt       = new TH2F("h_IEtaIPhiMapEvt", "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
  h_IEtaIPhiMapEvt->Sumw2();
  h_IEtaIPhiMapEvt_Ecut  = new TH2F("h_IEtaIPhiMapEvt_Ecut", "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
  h_IEtaIPhiMapEvt_Ecut->Sumw2();
  for (int i = 0; i < 50; i++){
    h_IEtaIPhiMapEvt_Cl[i]  = new TH2F(Form("h_IEtaIPhiMapEvt_Cl%d",i), "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
    h_IEtaIPhiMapEvt_Cl[i]->Sumw2();
    h_IEtaIPhiMapEvt_MCPart[i]  = new TH2F(Form("h_IEtaIPhiMapEvt_MCPart%d",i), "",  towersx+1, -0.5, towersx+0.5, towersy+1, -0.5, towersy+0.5);
    h_IEtaIPhiMapEvt_MCPart[i]->Sumw2();
  }

  TH2F* h_CylindricalMapEvt                 = new TH2F("h_CylindricalMapEvt", "", towersy, -TMath::Pi(), TMath::Pi(), towersx, zBins);
  h_CylindricalMapEvt->Sumw2();
  TH2F* h_CylindricalMapEvt_Cl[50]         = {NULL};      
  for(int i=0; i < 50; i++){
    h_CylindricalMapEvt_Cl[i]           = new TH2F(Form("h_CylindricalMapEvt_Cl%d",i), "", towersy, -TMath::Pi(), TMath::Pi(), towersx, zBins);
    h_CylindricalMapEvt_Cl[i]->Sumw2();
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
      case 6:           // HCAL in
        _nTowers    = _nTowers_HCALIN;
        _tower_E    = _tower_HCALIN_E;
        _tower_iEta = _tower_HCALIN_iEta;
        _tower_iPhi = _tower_HCALIN_iPhi;
        _tower_trueID = _tower_HCALIN_trueID;
        break;
      case 7:           // HCAL out
        _nTowers    = _nTowers_HCALOUT;
        _tower_E    = _tower_HCALOUT_E;
        _tower_iEta = _tower_HCALOUT_iEta;
        _tower_iPhi = _tower_HCALOUT_iPhi;
        _tower_trueID = _tower_HCALOUT_trueID;
        break;
      case 10:          // BECAL 
        _nTowers    = _nTowers_BECAL;
        _tower_E    = _tower_BECAL_E;
        _tower_iEta = _tower_BECAL_iEta;
        _tower_iPhi = _tower_BECAL_iPhi;
        _tower_trueID = _tower_BECAL_trueID;
        break;
      default:
        std::cout << "Wrong caloID chosen" << std::endl;
        return;
      }

    if(_nTowers == 0 ) return;

    for(int itwr=0;(int)itwr<_nTowers;itwr++){
      sumedE += _tower_E[itwr];
      if (_tower_E[itwr] > 1e-3 && verbosity) cout << _tower_iEta[itwr] << "\t" << _tower_iPhi[itwr] << "\t"<< _tower_E[itwr] << endl;
      h_IEtaIPhiMapEvt->Fill(_tower_iEta[itwr],_tower_iPhi[itwr],_tower_E[itwr]);
      if(_tower_E[itwr]>minECutOff){
        h_IEtaIPhiMapEvt_Ecut->Fill(_tower_iEta[itwr],_tower_iPhi[itwr],_tower_E[itwr]);
      }
      float tmpZ    = arr_GeoStrt[_tower_iEta[itwr]-1][_tower_iPhi[itwr]-1].tower_z;
      float tmpPhi  = arr_GeoStrt[_tower_iEta[itwr]-1][_tower_iPhi[itwr]-1].tower_Phi;
    //   std::cout << tmpZ << "\t" << tmpPhi << "\t" << _tower_E[itwr] << std::endl;
      h_CylindricalMapEvt->Fill(tmpPhi,tmpZ,_tower_E[itwr]);
    }
  }

  int nclusters = 0;

  // vector of towers used in the clusterizer
  std::vector<towersStrct> input_towers;
  std::vector<towersStrct> input_towersForMC;
  // vector of towers within the currently found cluster
  std::vector<towersStrct> cluster_towers;
  std::vector<clustersStrct> rec_clusters;
  
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
      if (verbosity) std::cout << tempstructT.tower_E << "\t" << tempstructT.tower_iEta << "\t" << tempstructT.tower_iPhi << std::endl;
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
        float eTwr    = cluster_towers.at(tic).tower_E;
        float tmpZ    = arr_GeoStrt[iEtaTwr-1][iPhiTwr-1].tower_z;
        float tmpPhi  = arr_GeoStrt[iEtaTwr-1][iPhiTwr-1].tower_Phi;        
        if (verbosity) cout << iEtaTwr << "\t"<< iPhiTwr << "\t" << eTwr << endl;
        h_CylindricalMapEvt_Cl[nclusters]->Fill(tmpPhi,tmpZ,eTwr);
        h_IEtaIPhiMapEvt_Cl[nclusters]->Fill(iEtaTwr, iPhiTwr, eTwr);
      }
      
      rec_clusters.push_back(tempstructC);

      nclusters++;
      
    } else {
      input_towers.clear();
    }
  }

  cout << "here" << endl;
  int nMCCurr   = 0;
  int nMCCurrID = input_towersForMC.at(0).tower_trueID;
  while (!input_towersForMC.empty() && nMCCurr < 49) {
    if (nMCCurrID != input_towersForMC.at(0).tower_trueID){
      nMCCurr++;
      nMCCurrID = input_towersForMC.at(0).tower_trueID;
    }

    int iEtaTwr = input_towersForMC.at(0).tower_iEta;
    int iPhiTwr = input_towersForMC.at(0).tower_iPhi;
    float eTwr    = input_towersForMC.at(0).tower_E;
    float tmpZ    = arr_GeoStrt[iEtaTwr-1][iPhiTwr-1].tower_z;
    float tmpPhi  = arr_GeoStrt[iEtaTwr-1][iPhiTwr-1].tower_Phi;        
    h_IEtaIPhiMapEvt_MCPart[nMCCurr]->Fill(iEtaTwr, iPhiTwr, eTwr);
    input_towersForMC.erase(input_towersForMC.begin());
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
  Int_t pixelBaseCanvas = 1000;
  Float_t fracX = (Float_t)towersx/maxTower;
  Float_t fracY = (Float_t)towersy/maxTower;
  Int_t pixelX, pixelY;
  Float_t marginBottom, marginTop, marginLeft, marginRight;
  Float_t labelOffset, tickLength;

  switch( caloID ){
      case 6:           // HCAL in
        pixelX  = 35*3.4+fracX*pixelBaseCanvas;
        pixelY  = 35*1.6+fracY*pixelBaseCanvas;
        marginBottom = (Float_t)(35*1.5)/(Float_t)pixelY;
        marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
        marginLeft = 35*1.5/(Float_t)pixelX;
        marginRight = 35*1.9/(Float_t)pixelX;
        labelOffset = -0.015;
        tickLength = 0.02;
        break;
      case 7:           // HCAL out
        pixelX  = 35*3.4+fracX*pixelBaseCanvas;
        pixelY  = 35*1.6+fracY*pixelBaseCanvas;
        marginBottom = (Float_t)(35*1.5)/(Float_t)pixelY;
        marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
        marginLeft = 35*1.5/(Float_t)pixelX;
        marginRight = 35*1.9/(Float_t)pixelX;
        labelOffset = -0.015;
        tickLength = 0.02;
        break;
      case 10:          // BECAL 
        pixelX  = 35*4.9+fracX*pixelBaseCanvas;
        pixelY  = 35*2+fracY*pixelBaseCanvas;
        marginBottom = (Float_t)(35*1.9)/(Float_t)pixelY;
        marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
        marginLeft = 35*2.2/(Float_t)pixelX;
        marginRight = 35*2.7/(Float_t)pixelX;
        labelOffset = -0.01;
        tickLength = 0.02;
        break;
      default:
        std::cout << "Wrong caloID chosen" << std::endl;
        return;
  }
  
  cout << pixelX << "\t" << pixelY << "\t" << marginLeft << "\t" << marginRight << "\t" << marginTop << "\t" << marginBottom << endl;
  
  TCanvas* c2DColz = new TCanvas("c2DColz","",0,0,pixelX,pixelY);
  DrawGammaCanvasSettings( c2DColz, marginLeft, marginRight, marginTop, marginBottom);
  c2DColz->SetLogz(1);

  // 1D PLOT
  if( caloID == 10 ) SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_Ecut, "tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.58,1.05);
  else SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_Ecut, "tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.45,1.05);
  
  h_IEtaIPhiMapEvt_Ecut->GetXaxis()->SetLabelOffset(labelOffset);
  h_IEtaIPhiMapEvt_Ecut->GetXaxis()->SetTickLength(tickLength);
  
  h_IEtaIPhiMapEvt_Ecut->Draw("colz");
    
  drawLatexAdd(Form("%s, %s",labelDet.Data(), labelEnergy.Data()),1-marginRight-0.05,1-marginTop-0.05,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(labelParticlProp,1-marginRight-0.05,1-marginTop-0.05-0.9*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  if (collisionsSys.Contains("Single")){
    drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),1-marginRight-0.05,1-marginTop-0.05-2*0.9*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  }

  c2DColz->Print(Form("%s/2D/2D_Event%d_Ecut.%s", outputDirLocal.Data(),plotEvent, suffix.Data()));

  // reset 2D canvas


  switch( caloID ){
      case 6:           // HCAL in
        pixelX  = 35*2.7+fracX*pixelBaseCanvas;
        pixelY  = 35*1.5+fracY*pixelBaseCanvas;
        marginBottom = (Float_t)(35*1.4)/(Float_t)pixelY;
        marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
        marginLeft = 35*1.4/(Float_t)pixelX;
        marginRight = 35*0.3/(Float_t)pixelX;
        labelOffset = -0.015;
        tickLength = 0.02;
        break;
      case 7:           // HCAL out
        pixelX  = 35*2.7+fracX*pixelBaseCanvas;
        pixelY  = 35*1.5+fracY*pixelBaseCanvas;
        marginBottom = (Float_t)(35*1.4)/(Float_t)pixelY;
        marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
        marginLeft = 35*1.4/(Float_t)pixelX;
        marginRight = 35*0.3/(Float_t)pixelX;
        labelOffset = -0.015;
        tickLength = 0.02;
        break;
      case 10:          // BECAL 
        pixelX  = 35*2.2+fracX*pixelBaseCanvas;
        pixelY  = 35*2+fracY*pixelBaseCanvas;
        marginBottom = (Float_t)(35*1.9)/(Float_t)pixelY;
        marginTop = (Float_t)(35*0.1)/(Float_t)pixelY;
        marginLeft = 35*1.9/(Float_t)pixelX;
        marginRight = 35*0.3/(Float_t)pixelX;
        labelOffset = -0.01;
        tickLength = 0.02;
        break;
      default:
        std::cout << "Wrong caloID chosen" << std::endl;
        return;
  }

  TCanvas* c2DBox = new TCanvas("c2DBox","",0,0,pixelX,pixelY);

  DrawGammaCanvasSettings( c2DBox, marginLeft, marginRight, marginTop, marginBottom);
  
  //************************************************************************************
  // draw 2D hist with box option
  //************************************************************************************  
  if( caloID != 10 )   SetStyleHistoTH2ForGraphs(h_IEtaIPhiMapEvt_Ecut, "tower #eta index", "tower #phi index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.4,1.0);

  h_IEtaIPhiMapEvt_Ecut->SetLineColor(kGray+2);
  h_IEtaIPhiMapEvt_Ecut->Draw("box");
  
  drawLatexAdd(Form("%s, %s",labelDet.Data(), labelEnergy.Data()),1-marginRight-0.05,marginBottom+0.03+0.9*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
  drawLatexAdd(labelParticlProp,1-marginRight-0.05,marginBottom+0.03,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
  if (collisionsSys.Contains("Single")){
    drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),marginLeft+0.04,1-marginTop-textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  }
  c2DBox->Print(Form("%s/2D/2D_Event%d_Ecut_Box.%s", outputDirLocal.Data(),plotEvent, suffix.Data()));
  
  //************************************************************************************
  // draw 2D hist with box option and different clusters on top
  //************************************************************************************
  for (int nCl = 0; nCl < nclusters; nCl++){
    h_IEtaIPhiMapEvt_Cl[nCl]->SetLineColor(colorCluster[nCl]);
    h_IEtaIPhiMapEvt_Cl[nCl]->Draw("box,same");
  }
  if (collisionsSys.Contains("Single")){
    if(nclusters>0)drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", firstClusterE ),marginLeft+0.04,1-marginTop-2*0.85*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  }
  c2DBox->Print(Form("%s/2DwClusters/2D_Event%d_Ecut_clusters.%s", outputDirLocal.Data(),plotEvent, suffix.Data()));

  //************************************************************************************
  // draw 2D hist with box option and different MC particles on top
  //************************************************************************************  
  h_IEtaIPhiMapEvt_Ecut->SetLineColor(kGray+2);
  h_IEtaIPhiMapEvt_Ecut->Draw("box");
  
  for (int nMC = 0; nMC < nMCCurr; nMC++){
    h_IEtaIPhiMapEvt_MCPart[nMC]->SetLineColor(colorCluster[nMC]);
    h_IEtaIPhiMapEvt_MCPart[nMC]->Draw("box,same");
  }
  if (collisionsSys.Contains("Single")){
    if(nclusters>0)drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", firstClusterE ),marginLeft+0.04,1-marginTop-2*0.85*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  }
  drawLatexAdd(Form("%s, %s",labelDet.Data(), labelEnergy.Data()),1-marginRight-0.05,marginBottom+0.03+0.9*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
  drawLatexAdd(labelParticlProp,1-marginRight-0.05,marginBottom+0.03,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);  
  c2DBox->Print(Form("%s/2DwClusters/2D_Event%d_Ecut_MCPart.%s", outputDirLocal.Data(),plotEvent, suffix.Data()));

  //************************************************************************************
  // draw 3D plot in cylindrical coordinates
  //************************************************************************************

  TCanvas* c3DCyl = new TCanvas("c3DCyl","",0,0,400,300);
  DrawGammaCanvasSettings( c3DCyl, 0, 0, 0, 0);
  c3DCyl->SetTheta(60);

  for(int nCl=0; nCl < nclusters; nCl++){
    h_CylindricalMapEvt->Add(h_CylindricalMapEvt_Cl[nCl],-1);
  }

  THStack *h_CylStack    = new THStack("h_CylStack","");
  h_CylindricalMapEvt->SetLineColor(kGray);
  h_CylStack->Add(h_CylindricalMapEvt);

  if(nclusters>1) std::cout << "Event: " << plotEvent << std::endl;
  for(int nCl=0; nCl < nclusters; nCl++){
    h_CylindricalMapEvt_Cl[nCl]->SetLineColor(colorCluster[nCl]);
    h_CylindricalMapEvt_Cl[nCl]->SetMaximum(h_CylindricalMapEvt->GetMaximum());
    h_CylStack->Add(h_CylindricalMapEvt_Cl[nCl]);
  }
  h_CylStack->Draw("LEGO1 CYL");


  drawLatexAdd(Form("%s, %s",labelDet.Data(), labelEnergy.Data()),0.03,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  drawLatexAdd(labelParticlProp,0.03,0.92-0.9*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  if (collisionsSys.Contains("Single")){
    drawLatexAdd(Form("#it{E}_{calo} = %0.1f GeV",sumedE ),0.03,0.92-2*0.9*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    if(nclusters>0)drawLatexAdd(Form("#it{E}_{cl} = %0.1f GeV", firstClusterE ),0.8,0.92-2*0.85*textSizeLabelsRel,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  }

  c3DCyl->Print(Form("%s/3D/CylTest_Lego_%d.pdf",outputDirLocal.Data(),plotEvent));
  
}
