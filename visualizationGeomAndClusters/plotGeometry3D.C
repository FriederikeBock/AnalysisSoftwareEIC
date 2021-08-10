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


void plotGeometry3D(
    TString nameFile          = "",
    TString suffix            = "pdf",
    TString detector          = "",
    Int_t caloIndex           = -1,
    bool do3D                 = false,
    Int_t debugLevel          = 0
    
){
  // dimensions 
  int towersx       = 600;
  int towersy       = 600;
  int towersz       = 1050;

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  
  TString outputDir                             = Form("plots/%s",detector.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  double textSizeSinglePad                    = 0.05;
  double textSizeLabelsPixel                  = 35;
  double textSizeLabelsRel                    = 58./1300;

  int granularity = 1;
  TH3F* h_XYZ                               = new TH3F("h_XYZ", "", towersz, -500, 550,  towersx, -300, 300, towersy, -300, 300);
  TH3F* h_XYZ_iEta                          = new TH3F("h_XYZ_iEta", "", towersz, -500, 550,  towersx, -300, 300, towersy, -300, 300);
  TH3F* h_XYZ_iPhi                          = new TH3F("h_XYZ_IPhi", "", towersz, -500, 550,  towersx, -300, 300, towersy, -300, 300);
  TH2F* h_PhiEta                            = new TH2F("h_EtaPhi", "", 800, -4, 4, 400, -TMath::Pi(), TMath::Pi());
  TH2F* h_PhiEta_iEta                       = new TH2F("h_PhiEta_IEta", "", 800, -4, 4, 400, -TMath::Pi(), TMath::Pi());
  TH2F* h_PhiEta_iPhi                       = new TH2F("h_PhiEta_IPhi", "", 800, -4, 4, 400, -TMath::Pi(), TMath::Pi());
  TH2F* h_RZ                                = new TH2F("h_ZR", "", towersz, -500, 550, 600, 0, 300);
  TH2F* h_ZR_iEta                           = new TH2F("h_ZR_iEta", "", towersz, -500, 550, 600, 0, 300);
  TH2F* h_ZR_iPhi                           = new TH2F("h_ZR_iPhi", "", towersz, -500, 550, 600, 0, 300);
  
  // read file
  
  const int _maxNTowersCalo = 5000000;
  int _calogeom_ID;
  int _calogeom_towers_N;
  int*  _calogeom_towers_iEta = new int[_maxNTowersCalo];
  int*  _calogeom_towers_iPhi = new int[_maxNTowersCalo];
  int*  _calogeom_towers_iL   = new int[_maxNTowersCalo];
  float*  _calogeom_towers_Eta = new float[_maxNTowersCalo];
  float*  _calogeom_towers_Phi = new float[_maxNTowersCalo];
  float*  _calogeom_towers_x = new float[_maxNTowersCalo];
  float*  _calogeom_towers_y = new float[_maxNTowersCalo];
  float*  _calogeom_towers_z = new float[_maxNTowersCalo];

  float minX = 10000;
  float minY = 10000;
  float minZ = 10000;
  float maxX = -10000;
  float maxY = -10000;
  float maxZ = -10000;
  float minIEta = 1000;
  float maxIEta = -1000;
  float minIPhi = 1000;
  float maxIPhi = -1000;
  float minIZ = 1000;
  float maxIZ = -1000;
  
  TChain * tt_geometry = new TChain("geometry_tree", "geometry_tree");
  if(tt_geometry){
    tt_geometry->Add(nameFile.Data());
    tt_geometry->SetBranchAddress("calo",              &_calogeom_ID);
    tt_geometry->SetBranchAddress("calo_towers_N",     &_calogeom_towers_N);
    tt_geometry->SetBranchAddress("calo_towers_iEta",  _calogeom_towers_iEta);
    tt_geometry->SetBranchAddress("calo_towers_iPhi",  _calogeom_towers_iPhi);
    tt_geometry->SetBranchAddress("calo_towers_iL",    _calogeom_towers_iL);
    tt_geometry->SetBranchAddress("calo_towers_Eta",   _calogeom_towers_Eta);
    tt_geometry->SetBranchAddress("calo_towers_Phi",   _calogeom_towers_Phi);
    tt_geometry->SetBranchAddress("calo_towers_x",     _calogeom_towers_x);
    tt_geometry->SetBranchAddress("calo_towers_y",     _calogeom_towers_y);
    tt_geometry->SetBranchAddress("calo_towers_z",     _calogeom_towers_z);
    
    Int_t nentries = Int_t(tt_geometry->GetEntries());
    cout << nentries << endl;
    for (Long64_t i=0; i<nentries;i++) {
      tt_geometry->GetEntry(i);
      minIEta = 1000;
      maxIEta = -1000;
      minIPhi = 1000;
      maxIPhi = -1000;
      minIZ = 1000;
      maxIZ = -1000;
    
      if (_calogeom_ID == caloIndex && caloIndex != -1){
        cout << "found calo: "<< _calogeom_ID << "\t"<< _calogeom_towers_N << endl;
        for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
          if (_calogeom_towers_x[itow] < minX)
            minX = _calogeom_towers_x[itow];
          if (_calogeom_towers_y[itow] < minY)
            minY = _calogeom_towers_y[itow];
          if (_calogeom_towers_z[itow] < minZ)
            minZ = _calogeom_towers_z[itow];
          if (_calogeom_towers_x[itow] > maxX)
            maxX = _calogeom_towers_x[itow];
          if (_calogeom_towers_y[itow] > maxY)
            maxY = _calogeom_towers_y[itow];
          if (_calogeom_towers_z[itow] > maxZ)
            maxZ = _calogeom_towers_z[itow];
          
          if (_calogeom_towers_iEta[itow] < minIEta)
            minIEta = _calogeom_towers_iEta[itow];
          if (_calogeom_towers_iPhi[itow] < minIPhi)
            minIPhi = _calogeom_towers_iPhi[itow];
          if (_calogeom_towers_iL[itow] < minIZ)
            minIZ = _calogeom_towers_iL[itow];
          
          if (_calogeom_towers_iEta[itow] > maxIEta)
            maxIEta = _calogeom_towers_iEta[itow];
          if (_calogeom_towers_iPhi[itow] > maxIPhi)
            maxIPhi = _calogeom_towers_iPhi[itow];
          if (_calogeom_towers_iL[itow] > maxIZ)
            maxIZ = _calogeom_towers_iL[itow];
          
          h_XYZ->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  1);
          h_PhiEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 1);
          h_RZ->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 1);
          if (_calogeom_towers_iEta[itow] == 30 || _calogeom_towers_iEta[itow] == 70 || _calogeom_towers_iEta[itow] == 1){
            h_XYZ_iEta->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  20);
            h_PhiEta_iEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 20);
            h_ZR_iEta->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 20);
            if(debugLevel) cout << itow << "\t" <<  _calogeom_towers_x[itow] << "\t" <<  _calogeom_towers_y[itow] << "\t" <<  _calogeom_towers_z[itow] << "\t" <<  _calogeom_towers_Eta[itow] << "\t" <<  _calogeom_towers_Phi[itow]  << endl;
          }
          if (_calogeom_towers_iPhi[itow] == 1 || _calogeom_towers_iPhi[itow] == 50 || _calogeom_towers_iPhi[itow] == 100){
            h_XYZ_iPhi->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  10);
            h_PhiEta_iPhi->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 10);            
            h_ZR_iPhi->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 10);
          }
        }
      } else if (caloIndex == -1){
        if (_calogeom_ID == 1 || _calogeom_ID == 3 || _calogeom_ID == 4|| _calogeom_ID == 9 || _calogeom_ID == 10){
          cout << "found calo: "<< _calogeom_ID << "\t"<< _calogeom_towers_N << endl;
          for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
            if (_calogeom_towers_x[itow] < minX)
              minX = _calogeom_towers_x[itow];
            if (_calogeom_towers_y[itow] < minY)
              minY = _calogeom_towers_y[itow];
            if (_calogeom_towers_z[itow] < minZ)
              minZ = _calogeom_towers_z[itow];
            if (_calogeom_towers_x[itow] > maxX)
              maxX = _calogeom_towers_x[itow];
            if (_calogeom_towers_y[itow] > maxY)
              maxY = _calogeom_towers_y[itow];
            if (_calogeom_towers_z[itow] > maxZ)
              maxZ = _calogeom_towers_z[itow];
            
            if (_calogeom_towers_iEta[itow] < minIEta)
              minIEta = _calogeom_towers_iEta[itow];
            if (_calogeom_towers_iPhi[itow] < minIPhi)
              minIPhi = _calogeom_towers_iPhi[itow];
            if (_calogeom_towers_iL[itow] < minIZ)
              minIZ = _calogeom_towers_iL[itow];
            
            if (_calogeom_towers_iEta[itow] > maxIEta)
              maxIEta = _calogeom_towers_iEta[itow];
            if (_calogeom_towers_iPhi[itow] > maxIPhi)
              maxIPhi = _calogeom_towers_iPhi[itow];
            if (_calogeom_towers_iL[itow] > maxIZ)
              maxIZ = _calogeom_towers_iL[itow];
            
            h_XYZ->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  1);
            h_PhiEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 1);
            h_RZ->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 1);
            if (_calogeom_towers_iEta[itow] == 30 || _calogeom_towers_iEta[itow] == 70 || _calogeom_towers_iEta[itow] == 1){
              h_XYZ_iEta->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  20);
              h_PhiEta_iEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 20);
              h_ZR_iEta->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 20);
              if(debugLevel) cout << itow << "\t" <<  _calogeom_towers_x[itow] << "\t" <<  _calogeom_towers_y[itow] << "\t" <<  _calogeom_towers_z[itow] << "\t" <<  _calogeom_towers_Eta[itow] << "\t" <<  _calogeom_towers_Phi[itow]  << endl;
            }
            if (_calogeom_towers_iPhi[itow] == 1 || _calogeom_towers_iPhi[itow] == 50 || _calogeom_towers_iPhi[itow] == 100){
              h_XYZ_iPhi->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  10);
              h_PhiEta_iPhi->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 10);            
              h_ZR_iPhi->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 10);
            }
          }
        } else {
          continue;
        }
      } else if (caloIndex == -2){
        if (_calogeom_ID == 0 || _calogeom_ID == 2 || _calogeom_ID == 5|| _calogeom_ID == 6 || _calogeom_ID == 7 || _calogeom_ID == 8){
          cout << "found calo: "<< _calogeom_ID << "\t"<< _calogeom_towers_N << endl;
          for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
            if (_calogeom_towers_x[itow] < minX)
              minX = _calogeom_towers_x[itow];
            if (_calogeom_towers_y[itow] < minY)
              minY = _calogeom_towers_y[itow];
            if (_calogeom_towers_z[itow] < minZ)
              minZ = _calogeom_towers_z[itow];
            if (_calogeom_towers_x[itow] > maxX)
              maxX = _calogeom_towers_x[itow];
            if (_calogeom_towers_y[itow] > maxY)
              maxY = _calogeom_towers_y[itow];
            if (_calogeom_towers_z[itow] > maxZ)
              maxZ = _calogeom_towers_z[itow];

            if (_calogeom_towers_iEta[itow] < minIEta)
              minIEta = _calogeom_towers_iEta[itow];
            if (_calogeom_towers_iPhi[itow] < minIPhi)
              minIPhi = _calogeom_towers_iPhi[itow];
            if (_calogeom_towers_iL[itow] < minIZ)
              minIZ = _calogeom_towers_iL[itow];
            
            if (_calogeom_towers_iEta[itow] > maxIEta)
              maxIEta = _calogeom_towers_iEta[itow];
            if (_calogeom_towers_iPhi[itow] > maxIPhi)
              maxIPhi = _calogeom_towers_iPhi[itow];
            if (_calogeom_towers_iL[itow] > maxIZ)
              maxIZ = _calogeom_towers_iL[itow];
            
            h_XYZ->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  1);
            h_PhiEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 1);
            h_RZ->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 1);
            if (_calogeom_towers_iEta[itow] == 30 || _calogeom_towers_iEta[itow] == 70 || _calogeom_towers_iEta[itow] == 1){
              h_XYZ_iEta->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  20);
              h_PhiEta_iEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 20);
              h_ZR_iEta->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 20);
              if(debugLevel) cout << itow << "\t" <<  _calogeom_towers_x[itow] << "\t" <<  _calogeom_towers_y[itow] << "\t" <<  _calogeom_towers_z[itow] << "\t" <<  _calogeom_towers_Eta[itow] << "\t" <<  _calogeom_towers_Phi[itow]  << endl;
            }
            if (_calogeom_towers_iPhi[itow] == 1 || _calogeom_towers_iPhi[itow] == 50 || _calogeom_towers_iPhi[itow] == 100){
              h_XYZ_iPhi->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  10);
              h_PhiEta_iPhi->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 10);            
              h_ZR_iPhi->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 10);
            }
          }
        } else {
          continue;
        }
      } else if (caloIndex == -3){
        cout << "found calo: "<< _calogeom_ID << "\t"<< _calogeom_towers_N << endl;
        for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
          if (_calogeom_towers_x[itow] < minX)
            minX = _calogeom_towers_x[itow];
          if (_calogeom_towers_y[itow] < minY)
            minY = _calogeom_towers_y[itow];
          if (_calogeom_towers_z[itow] < minZ)
            minZ = _calogeom_towers_z[itow];
          if (_calogeom_towers_x[itow] > maxX)
            maxX = _calogeom_towers_x[itow];
          if (_calogeom_towers_y[itow] > maxY)
            maxY = _calogeom_towers_y[itow];
          if (_calogeom_towers_z[itow] > maxZ)
            maxZ = _calogeom_towers_z[itow];

          if (_calogeom_towers_iEta[itow] < minIEta)
            minIEta = _calogeom_towers_iEta[itow];
          if (_calogeom_towers_iPhi[itow] < minIPhi)
            minIPhi = _calogeom_towers_iPhi[itow];
          if (_calogeom_towers_iL[itow] < minIZ)
            minIZ = _calogeom_towers_iL[itow];
          
          if (_calogeom_towers_iEta[itow] > maxIEta)
            maxIEta = _calogeom_towers_iEta[itow];
          if (_calogeom_towers_iPhi[itow] > maxIPhi)
            maxIPhi = _calogeom_towers_iPhi[itow];
          if (_calogeom_towers_iL[itow] > maxIZ)
            maxIZ = _calogeom_towers_iL[itow];

          h_XYZ->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  1);
          h_PhiEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 1);
          h_RZ->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 1);
          if (_calogeom_towers_iEta[itow] == 30 || _calogeom_towers_iEta[itow] == 70 || _calogeom_towers_iEta[itow] == 1){
            h_XYZ_iEta->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  20);
            h_PhiEta_iEta->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 20);
            h_ZR_iEta->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 20);
            if(debugLevel) cout << itow << "\t" <<  _calogeom_towers_x[itow] << "\t" <<  _calogeom_towers_y[itow] << "\t" <<  _calogeom_towers_z[itow] << "\t" <<  _calogeom_towers_Eta[itow] << "\t" <<  _calogeom_towers_Phi[itow]  << endl;
          }
          if (_calogeom_towers_iPhi[itow] == 1 || _calogeom_towers_iPhi[itow] == 50 || _calogeom_towers_iPhi[itow] == 100){
            h_XYZ_iPhi->Fill(_calogeom_towers_z[itow], _calogeom_towers_x[itow], _calogeom_towers_y[itow],  10);
            h_PhiEta_iPhi->Fill( _calogeom_towers_Eta[itow], _calogeom_towers_Phi[itow], 10);            
            h_ZR_iPhi->Fill(  _calogeom_towers_z[itow], TMath::Sqrt(_calogeom_towers_x[itow]*_calogeom_towers_x[itow]+_calogeom_towers_y[itow]*_calogeom_towers_y[itow]), 10);
          }
        }
      
      } else {
        continue;
      }
      cout << "iPhi: "<< minIPhi << "\t"<< maxIPhi << endl;
      cout << "iEta: "<< minIEta << "\t"<< maxIEta << endl;
      cout << "iL: "<< minIZ << "\t"<< maxIZ << endl;
    }
  }
  
  TCanvas* cReso = new TCanvas("cReso","",0,0,1000,950);
  if (do3D){
    //************************************************************************************
    // draw 3D hist with box option   
    //************************************************************************************
    DrawGammaCanvasSettings( cReso, 0.12, 0.05, 0.01, 0.1);
    cReso->SetTheta(18);
    cReso->SetPhi(40);
    
      SetStyleHistoTH3ForGraphs(h_XYZ, "z (cm)", "x (cm)", "y (cm)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1, 1.1, 1.15, 505, 510,510);

      h_XYZ->SetMaximum(20);
      h_XYZ->SetLineColor(kGray);
      h_XYZ->SetFillColor(kGray);
      h_XYZ->Draw("box");
      h_XYZ_iEta->SetLineColor(kBlue+2);
      h_XYZ_iEta->SetFillColor(kBlue+2);
      h_XYZ_iEta->Draw("box,same");
      h_XYZ_iPhi->SetLineColor(kRed+2);
      h_XYZ_iPhi->SetFillColor(kRed+2);
      h_XYZ_iPhi->Draw("box,same");

      drawLatexAdd(Form("ECCE %s", detector.Data()),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      
      TLegend* legendXYZ    = GetAndSetLegend2(0.01, 0.01, 0.4, 0.01+(3*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
      legendXYZ->AddEntry(h_XYZ, "all towers","lf");
      legendXYZ->AddEntry(h_XYZ_iEta, "towers iEta = 1, 30, 70","fl");
      legendXYZ->AddEntry(h_XYZ_iPhi, "towers iPhi = 1, 50, 100","fl");
      legendXYZ->Draw();
      
    cReso->Print(Form("%s/%s.%s", outputDir.Data(),detector.Data(), suffix.Data()));
      h_XYZ->SetMaximum(20);
      h_XYZ->GetXaxis()->SetRangeUser(minZ-10, maxZ+10);
      h_XYZ->GetYaxis()->SetRangeUser(minX-10, maxX+10);
      h_XYZ->GetZaxis()->SetRangeUser(minY-10, maxY+10);
      h_XYZ->Draw("box");
      h_XYZ_iEta->GetXaxis()->SetRangeUser(minZ-10, maxZ+10);
      h_XYZ_iEta->GetYaxis()->SetRangeUser(minX-10, maxX+10);
      h_XYZ_iEta->GetZaxis()->SetRangeUser(minY-10, maxY+10);
      h_XYZ_iEta->Draw("box,same");
      h_XYZ_iPhi->GetXaxis()->SetRangeUser(minZ-10, maxZ+10);
      h_XYZ_iPhi->GetYaxis()->SetRangeUser(minX-10, maxX+10);
      h_XYZ_iPhi->GetZaxis()->SetRangeUser(minY-10, maxY+10);
      h_XYZ_iPhi->Draw("box,same");

      drawLatexAdd(Form("ECCE %s", detector.Data()),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      
      legendXYZ->Draw();
      
    cReso->Print(Form("%s/%s_zoomed.%s", outputDir.Data(),detector.Data(), suffix.Data()));
  }
  
  DrawGammaCanvasSettings( cReso, 0.08, 0.02, 0.01, 0.08);
    SetStyleHistoTH2ForGraphs(h_PhiEta, "#eta", "#varphi", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.75, 0.73, 510,510);

    h_PhiEta->SetMaximum(20);
    h_PhiEta->SetLineColor(kBlack);
    h_PhiEta->SetFillColor(kBlack);
    h_PhiEta->Draw("box");
    h_PhiEta_iEta->SetLineColor(kBlue+2);
    h_PhiEta_iEta->SetFillColor(kBlue+2);
    h_PhiEta_iEta->Draw("box,same");
    h_PhiEta_iPhi->SetLineColor(kRed+2);
    h_PhiEta_iPhi->SetFillColor(kRed+2);
    h_PhiEta_iPhi->Draw("box,same");


    TLegend* legendPhiEta    = GetAndSetLegend2(0.12, 0.12, 0.4, 0.12+(3*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendPhiEta->AddEntry(h_PhiEta, "all towers","plf");
    legendPhiEta->AddEntry(h_PhiEta_iEta, "towers iEta = 1, 30, 70","fl");
    legendPhiEta->AddEntry(h_PhiEta_iPhi, "towers iPhi = 1, 50, 100","fl");
    legendPhiEta->Draw();
    
    drawLatexAdd(Form("ECCE %s", detector.Data()),0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
  cReso->Print(Form("%s/PhiEta%s.%s", outputDir.Data(),detector.Data(), suffix.Data()));

  DrawGammaCanvasSettings( cReso, 0.105, 0.01, 0.02, 0.08);
    SetStyleHistoTH2ForGraphs(h_RZ,  "Z (cm)", "R (cm)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.8, 1, 505,510);

    h_RZ->SetMaximum(20);
    h_RZ->SetLineColor(kBlack);
    h_RZ->SetFillColor(kBlack);
    h_RZ->Draw("box");
    h_ZR_iEta->SetLineColor(kBlue+2);
    h_ZR_iEta->SetFillColor(kBlue+2);
    h_ZR_iEta->Draw("box,same");
    h_ZR_iPhi->SetLineColor(kRed+2);
    h_ZR_iPhi->SetFillColor(kRed+2);
    h_ZR_iPhi->Draw("box,same");


    TLegend* legendRZ    = GetAndSetLegend2(0.12, 0.12, 0.4, 0.12+(3*0.85*textSizeLabelsRel),0.85*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendRZ->AddEntry(h_RZ, "all towers","plf");
    legendRZ->AddEntry(h_ZR_iEta, "towers iEta = 1, 30, 70","fl");
    legendRZ->AddEntry(h_ZR_iPhi, "towers iPhi = 1, 50, 100","fl");
    legendRZ->Draw();
    
    drawLatexAdd(Form("ECCE %s", detector.Data()),0.95,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    
  cReso->Print(Form("%s/RZ%s.%s", outputDir.Data(),detector.Data(), suffix.Data()));
  
}

