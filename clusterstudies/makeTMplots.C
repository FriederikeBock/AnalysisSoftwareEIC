#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

enum calotype {
  kFHCAL          = 0,
  kFEMC           = 1,
  kDRCALO         = 2,
  kEEMC           = 3,
  kCEMC           = 4,
  kEHCAL          = 5,
  kHCALIN         = 6,
  kHCALOUT        = 7,
  kLFHCAL         = 8,
  kEEMCG          = 9,
  kBECAL          = 10,
  kFOCAL          = 11
};
const int maxcalo = 12;
int calogeomindex[maxcalo] = {0};

TString str_calorimeter[maxcalo]      = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC", "EHCAL", "HCALIN", "HCALOUT", "LFHCAL", "EEMCG", "BECAL", "FOCAL"};
TString str_calorimeterPlot[maxcalo]  = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC", "EHCAL", "IHCAL", "OHCAL", "LFHCAL", "EEMCG", "BEMC", "FOCAL"};
const int _active_calo = 12;
bool _curr_active_calo[maxcalo]   = {false,    true,    false,   true,  false,   false,     true,      true,     true,   false,    true,  false};

int _combCalo[maxcalo]           = {kFEMC,  -1,     -1,       -1,      -1,      kEEMC,  kBECAL,   kBECAL,     kFEMC,    -1,       -1,       -1};
int _combCalo2[maxcalo]          = {-1,     -1,     -1,       -1,      -1,      -1,     -1,       -1,         -1,       -1,       -1,       -1};

enum clusterizertype {
    kMA         = 0,
    kV3         = 1,
    kV1         = 2,
    k5x5        = 3,
    kC5         = 4,
    kC3         = 5,
    k3x3        = 6,
    kDummy      = 7
};
// TString str_clusterizer[7] = {"V1", "V3", "3x3", "5x5", "C3", "C5", "MA"};
TString str_clusterizer[7] = {"MA", "V3", "V1", "5x5", "C5", "C3", "3x3"};
const int _active_algo = 1;


bool IsHCALCalorimeter(int caloID){
  switch (caloID){
    case kFEMC: return false;
    case kEEMC: return false;
    case kHCALIN: return true;
    case kHCALOUT: return true;
    case kLFHCAL: return true;
    case kBECAL: return false;
    default:
      std::cout << "IsHCALCalorimeter: caloID " << caloID << " not defined, returning false" << std::endl;
      return false;
  }
  return false;
}
bool IsForwardCalorimeter(int caloID){
  switch (caloID){
    case kFEMC: return true;
    case kEEMC: return true;
    case kHCALIN: return false;
    case kHCALOUT: return false;
    case kLFHCAL: return true;
    case kBECAL: return false;
    default:
      cout << "IsForwardCalorimeter: caloID " << caloID << " not defined, returning false" << endl;
      return false;
  }
  return false;
}


int GetCaloDirection(int caloID){
  switch (caloID){
    case kFEMC: return 2;
    case kEEMC: return 0;
    case kHCALIN: return 1;
    case kHCALOUT: return 1;
    case kLFHCAL: return 2;
    case kBECAL: return 1;
    default:
      std::cout << "GetCaloDirection: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}
float ReturnTrackMatchingWindowForCalo(int caloID){
  switch (caloID){
    case kFEMC: return 5.5;
    case kEEMC: return 4.0;
    case kHCALIN: return 0.1;
    case kHCALOUT: return 0.1;
    case kLFHCAL: return 10;
    case kBECAL: return 0.2;
    default:
      cout << "ReturnTrackMatchingWindowForCalo: caloID " << caloID << " not defined, returning -1" << endl;
      return -1;
  }
  return -1;
}


// ***********************************************************************************
// matching windows for Calo's going inward in dPhi
// ***********************************************************************************
float ReturnPhiCaloMatching(int caloID, int caloID2){
  switch (caloID){
    case kDRCALO: return 0.1;
    case kFHCAL: return 0.1;
    case kFEMC: return 0.1;
    case kEHCAL: return 0.1;
    case kEEMC: return 0.1;
    case kHCALIN: return 0.1;
    case kHCALOUT: 
      switch (caloID2){
        case kBECAL:
          return 0.08;
        case kHCALIN:
          return 0.02;
        default:
          return 0.1;
      }
    case kCEMC: return 0.1;
    case kEEMCG: return 0.1;
    case kLFHCAL: return 0.07;
    case kBECAL: return 0.05;
    case kFOCAL: return 0.1;
    default:
      std::cout << "ReturnPhiCaloMatching: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}

// ***********************************************************************************
// matching windows for Calo's going inward in dEta
// ***********************************************************************************
float ReturnEtaCaloMatching(int caloID, int caloID2){
  switch (caloID){
    case kDRCALO: return 0.1;
    case kFHCAL: return 0.1;
    case kFEMC: return 0.1;
    case kEHCAL: return 0.1;
    case kEEMC: return 0.1;
    case kHCALIN: return 0.1;
    case kHCALOUT: 
      switch (caloID2){
        case kBECAL:
          return 0.05;
        case kHCALIN:
          return 0.02;
        default:
          return 0.1;
      }
    case kCEMC: return 0.1;
    case kEEMCG: return 0.1;
    case kLFHCAL: return 0.1;
    case kBECAL: return 0.05;
    case kFOCAL: return 0.1;
    default:
      std::cout << "ReturnEtaCaloMatching: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}

void PalColor(){
    static Int_t  colors[50];
    static Bool_t initialized = kFALSE;

    Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[5]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[5] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[5]  = { 0.51, 1., 0.12, 0.00, 0.00};

    if(!initialized){
        Int_t FI = TColor::CreateGradientColorTable(5, stops, red, green, blue, 50);
        for (int i=0; i<50; i++) colors[i] = FI+i;
        initialized = kTRUE;
        return;
    }
    gStyle->SetPalette(50,colors);
}
void makeTMplots(
    TString nameFile  = "definebelow.root",
    TString suffix            = "pdf",
    bool plotUncut = false
){
  if(plotUncut)
    nameFile  = "/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPION_TTLGEO7_Match_opencut/output_TMSTUD.root";
  else
    nameFile  = "/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPION_TTLGEO7_Match/output_TMSTUD.root";
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput                       = ReturnDateStringForOutput();
  TString perfLabel                 = "#it{#bf{ECCE}} simulation";

  PalColor();
  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFile = new TFile(nameFile.Data());

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;


  // Double_t partEtaCalo[nEta+1]    = { -4.0, -3.5, -3.0, -2.5, -2.0,    -1.7, -1.2, -0.4, 0.4, 0.9, 
  //                                      1.3, 1.8, 2.5, 3.0, 3.5,        4.0};
  Int_t minEtaBinCaloTMS[3]       = {1, 5,  11};
  Int_t maxEtaBinCaloTMS[3]       = {5, 10, 15};

  TH2F*  h_TMstudies_2D_clus[_active_calo] = {nullptr};                         // [calorimeter_enum]

  TH1F*  h_TMstudies_clusSpec[_active_calo] = {nullptr};                        // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Matched[_active_calo] = {nullptr};                // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_MatchedCalo[_active_calo] = {nullptr};            // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_MatchedTrackAndCalo[_active_calo] = {nullptr};    // [calorimeter_enum][algorithm_enum]

  TH1F*  h_TMstudies_clusSpec_Eta[_active_calo][nEta+1] = {{nullptr}};          // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_MatchedCalo_Eta[_active_calo][nEta+1] = {{nullptr}};    // [calorimeter_enum][algorithm_enum]
  TH1F*  h_TMstudies_clusSpec_MatchedTrackAndCalo_Eta[_active_calo][nEta+1] = {{nullptr}};    // [calorimeter_enum][algorithm_enum]
  TH1F*  h_TMstudies_clusSpec_Matched_Eta[_active_calo][nEta+1] = {{nullptr}};  // [calorimeter_enum]

  TH1F*  h_TMstudies_clusSpec_Cut[_active_calo] = {nullptr};                        // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Matched_Cut[_active_calo] = {nullptr};                // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_MatchedCalo_Cut[_active_calo] = {nullptr};            // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_MatchedTrackAndCalo_Cut[_active_calo] = {nullptr};    // [calorimeter_enum][algorithm_enum]


  TH2F*  h_TMstudies_2D_delta[_active_calo] = {nullptr};                        // [calorimeter_enum]
  TH2F*  h_TMstudies_dx_vsEta[_active_calo] = {nullptr};                        // [calorimeter_enum]
  TH2F*  h_TMstudies_dy_vsEta[_active_calo] = {nullptr};                        // [calorimeter_enum]
  TH2F*  h_TMstudies_dx_vsClsE[_active_calo][nEta+1] = {{nullptr}};             // [calorimeter_enum]
  TH2F*  h_TMstudies_dy_vsClsE[_active_calo][nEta+1] = {{nullptr}};             // [calorimeter_enum]

  TH2F*  h_TMstudiesCalo_2D_delta[_active_calo][_active_calo] = {{nullptr}};            // [calorimeter_enum]
  TH2F*  h_TMstudiesCalo_dEta_vsEta[_active_calo][_active_calo] = {{nullptr}};          // [calorimeter_enum]
  TH2F*  h_TMstudiesCalo_dPhi_vsEta[_active_calo][_active_calo] = {{nullptr}};          // [calorimeter_enum]
  TH2F*  h_TMstudiesCalo_dEta_vsClsE[_active_calo][_active_calo][nEta+1] = {{{nullptr}}}; // [calorimeter_enum]
  TH2F*  h_TMstudiesCalo_dPhi_vsClsE[_active_calo][_active_calo][nEta+1] = {{{nullptr}}}; // [calorimeter_enum]

  for(int icalo=0;icalo<_active_calo;icalo++){
    if(!_curr_active_calo[icalo]) continue;
    bool isFwd = IsForwardCalorimeter(icalo);
    Int_t caloDir           = GetCaloDirection(icalo);
    TString str2DhistoName = isFwd ? "x_y" : "eta_phi";

    h_TMstudies_2D_clus[icalo] 	= (TH2F*) inputFile->Get(Form("%s/h_TMstudies_%s_clus_%s_MA",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()));
    if(!h_TMstudies_2D_clus[icalo]) cout << Form("%s/h_TMstudies_%s_clus_%s_MA",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()) << " not found" << endl;
    else cout << Form("%s/h_TMstudies_%s_clus_%s_MA",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()) << " found" << endl;
  
    h_TMstudies_clusSpec[icalo]             = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_clusSpec_Matched[icalo]     = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_Matched_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));

    h_TMstudies_2D_delta[icalo]            = (TH2F*) inputFile->Get(Form("%s/h_TMstudies_2D_delta_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_dx_vsEta[icalo]            = (TH2F*) inputFile->Get(Form("%s/h_TMstudies_dx_vsEta_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_dy_vsEta[icalo]            = (TH2F*) inputFile->Get(Form("%s/h_TMstudies_dy_vsEta_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));

    
    // add histos for HCal cluster to ECal cluster matching
    if (IsHCALCalorimeter(icalo)){
      h_TMstudies_clusSpec_MatchedCalo[icalo] = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_MatchedCalo_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
      
      h_TMstudies_clusSpec_MatchedTrackAndCalo[icalo]            = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_MatchedTrackAndCalo_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
      if(!h_TMstudies_clusSpec_MatchedTrackAndCalo[icalo]) cout << "could not find " << Form("%s/h_TMstudies_clusSpec_MatchedTrackAndCalo_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()) << endl;

      if (_combCalo[icalo] != -1){
        if (_curr_active_calo[_combCalo[icalo]]){
          // h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]    = (TH2F*) inputFile2->Get(Form("%s/h_clusterizer_all_2D_deltaCalo_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
          //                                                                         str_calorimeter[_combCalo[icalo]].Data()));
          h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]    = (TH2F*) inputFile->Get(Form("%s/h_TMstudiesCalo_2D_delta_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
                                                                                  str_calorimeter[_combCalo[icalo]].Data()));
          h_TMstudiesCalo_dEta_vsEta[icalo][_combCalo[icalo]]  = (TH2F*) inputFile->Get(Form("%s/h_TMstudiesCalo_dEta_vsEta_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
                                                                                str_calorimeter[_combCalo[icalo]].Data()));
          h_TMstudiesCalo_dPhi_vsEta[icalo][_combCalo[icalo]]  = (TH2F*) inputFile->Get(Form("%s/h_TMstudiesCalo_dPhi_vsEta_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
                                                                                str_calorimeter[_combCalo[icalo]].Data()));
        }
      }
      
      if (_combCalo2[icalo] != -1){
        if (_curr_active_calo[_combCalo2[icalo]]){
          // h_TMstudiesCalo_2D_delta[icalo][_combCalo2[icalo]]   = (TH2F*) inputFile2->Get(Form("%s/h_clusterizer_all_2D_deltaCalo_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
          //                                                                         str_calorimeter[_combCalo2[icalo]].Data()));
          h_TMstudiesCalo_2D_delta[icalo][_combCalo2[icalo]]   = (TH2F*) inputFile->Get(Form("%s/h_TMstudiesCalo_2D_delta_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
                                                                                  str_calorimeter[_combCalo2[icalo]].Data()));
          h_TMstudiesCalo_dEta_vsEta[icalo][_combCalo2[icalo]] = (TH2F*) inputFile->Get(Form("%s/h_TMstudiesCalo_dEta_vsEta_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
                                                                                str_calorimeter[_combCalo2[icalo]].Data()));
          h_TMstudiesCalo_dPhi_vsEta[icalo][_combCalo2[icalo]] = (TH2F*) inputFile->Get(Form("%s/h_TMstudiesCalo_dPhi_vsEta_%s_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data(),
                                                                                str_calorimeter[_combCalo2[icalo]].Data()));
        }
      }          
    }
    
    for (Int_t et = minEtaBinCaloTMS[caloDir]; et<maxEtaBinCaloTMS[caloDir]; et++){
      Double_t etaMin = partEtaCalo[0];
      Double_t etaMax = partEtaCalo[nEta];
      if (et < nEta){
        etaMin = partEtaCalo[et];
        etaMax = partEtaCalo[et+1];
      }
      // reconstructed  particles
      h_TMstudies_clusSpec_Eta[icalo][et]          = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_Eta_%1.1f_%1.1f_%s_MA",str_calorimeter[icalo].Data(),etaMin,etaMax,str_calorimeter[icalo].Data()));
      if(!h_TMstudies_clusSpec_Cut[icalo])h_TMstudies_clusSpec_Cut[icalo] = (TH1F*) h_TMstudies_clusSpec_Eta[icalo][et]->Clone(Form("h_TMstudies_clusSpec_Cut_%s_MA",str_calorimeter[icalo].Data()));
      else h_TMstudies_clusSpec_Cut[icalo]->Add(h_TMstudies_clusSpec_Eta[icalo][et]);
      // reconstructed & matched particles

      h_TMstudies_clusSpec_Matched_Eta[icalo][et]  = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_Matched_Eta_%1.1f_%1.1f_%s_MA",str_calorimeter[icalo].Data(),etaMin,etaMax,str_calorimeter[icalo].Data()));
      if(!h_TMstudies_clusSpec_Matched_Cut[icalo])h_TMstudies_clusSpec_Matched_Cut[icalo] = (TH1F*) h_TMstudies_clusSpec_Matched_Eta[icalo][et]->Clone(Form("h_TMstudies_clusSpec_Matched_Cut_%s_MA",str_calorimeter[icalo].Data()));
      else h_TMstudies_clusSpec_Matched_Cut[icalo]->Add(h_TMstudies_clusSpec_Matched_Eta[icalo][et]);

      if (IsHCALCalorimeter(icalo)){
        h_TMstudies_clusSpec_MatchedCalo_Eta[icalo][et] = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_MatchedCalo_Eta_%1.1f_%1.1f_%s_MA",str_calorimeter[icalo].Data(),etaMin,etaMax,str_calorimeter[icalo].Data()));
        if(!h_TMstudies_clusSpec_MatchedCalo_Cut[icalo])h_TMstudies_clusSpec_MatchedCalo_Cut[icalo] = (TH1F*) h_TMstudies_clusSpec_MatchedCalo_Eta[icalo][et]->Clone(Form("h_TMstudies_clusSpec_MatchedCalo_Cut_%s_MA",str_calorimeter[icalo].Data()));
        else h_TMstudies_clusSpec_MatchedCalo_Cut[icalo]->Add(h_TMstudies_clusSpec_MatchedCalo_Eta[icalo][et]);

        h_TMstudies_clusSpec_MatchedTrackAndCalo_Eta[icalo][et] = (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_MatchedTrackAndCalo_Eta_%1.1f_%1.1f_%s_MA",str_calorimeter[icalo].Data(),etaMin,etaMax,str_calorimeter[icalo].Data()));
        if(!h_TMstudies_clusSpec_MatchedTrackAndCalo_Cut[icalo])h_TMstudies_clusSpec_MatchedTrackAndCalo_Cut[icalo] = (TH1F*) h_TMstudies_clusSpec_MatchedTrackAndCalo_Eta[icalo][et]->Clone(Form("h_TMstudies_clusSpec_MatchedTrackAndCalo_Cut_%s_MA",str_calorimeter[icalo].Data()));
        else h_TMstudies_clusSpec_MatchedTrackAndCalo_Cut[icalo]->Add(h_TMstudies_clusSpec_MatchedTrackAndCalo_Eta[icalo][et]);
      }
    }
  }
  Color_t linecolx = kGray;
  // 2D PLOT
  TCanvas* cXYplot = new TCanvas("cXYplot","",0,0,1000,1000);
  DrawGammaCanvasSettings( cXYplot, 0.13, 0.015, 0.015, 0.1);
  // cXYplot->SetLogz();
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(!_curr_active_calo[icalo]) continue;
    if(!h_TMstudies_2D_clus[icalo]) continue;
    cout << "plotting calo " << str_calorimeterPlot[icalo].Data() << endl;
    bool isFwd = IsForwardCalorimeter(icalo);
    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudies_2D_clus[icalo],"#it{#eta}", "#it{#varphi}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    // h_TMstudies_2D_clus[icalo]->GetXaxis()->SetRangeUser(-249,249);
    // h_TMstudies_2D_clus[icalo]->GetYaxis()->SetRangeUser(-249,249);
    h_TMstudies_2D_clus[icalo]->Draw("col");
    drawLatexAdd(Form("#pi^{#pm} clusters in %s",str_calorimeterPlot[icalo].Data()),0.17,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    cXYplot->Print(Form("%s/2D_clusters_%s.%s", outputDir.Data(),str_calorimeterPlot[icalo].Data(), suffix.Data()));
  }
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(!_curr_active_calo[icalo]) continue;
    if(!IsHCALCalorimeter(icalo)) continue;
    cout << "plotting calo " << str_calorimeterPlot[icalo].Data() << endl;
    bool isFwd = IsForwardCalorimeter(icalo);
    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]],Form("#it{#eta}_{%s}-#it{#eta}_{%s}",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()), Form("#it{#varphi}_{%s}-#it{#varphi}_{%s}",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()), 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]->GetXaxis()->SetRangeUser(-0.24,0.24);
    h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]->GetYaxis()->SetRangeUser(-0.24,0.24);
    h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]->Draw("col");
    float matchingwindowphi = ReturnPhiCaloMatching(icalo, _combCalo[icalo]);
    float matchingwindoweta = ReturnEtaCaloMatching(icalo, _combCalo[icalo]);
    if(plotUncut){
      DrawGammaLines(-matchingwindoweta, -matchingwindoweta, -matchingwindowphi, matchingwindowphi, 3, linecolx, 2);
      DrawGammaLines(matchingwindoweta, matchingwindoweta, -matchingwindowphi, matchingwindowphi,   3, linecolx, 2);
      DrawGammaLines(-matchingwindoweta, matchingwindoweta, -matchingwindowphi, -matchingwindowphi, 3, linecolx, 2);
      DrawGammaLines(-matchingwindoweta, matchingwindoweta, matchingwindowphi, matchingwindowphi,   3, linecolx, 2);
    }
    //  TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
    //  el4->SetFillStyle(0);
    //  el4->SetLineColor(kGray+2);
    //  el4->SetLineWidth(2);
    //  el4->SetLineStyle(7);
    //  el4->Draw("same");
      // drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%s and %s #pi^{#pm} clusters matched",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),0.17,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cXYplot->Print(Form("%s/2D_caloMatching_delta_%s_%s.%s", outputDir.Data(),str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data(), suffix.Data()));
  }
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(!_curr_active_calo[icalo]) continue;
    if(!IsHCALCalorimeter(icalo)) continue;
    cout << "plotting calo " << str_calorimeterPlot[icalo].Data() << endl;
    bool isFwd = IsForwardCalorimeter(icalo);
    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudiesCalo_dEta_vsEta[icalo][_combCalo[icalo]],Form("#it{#eta}_{%s}-#it{#eta}_{%s}",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()), "#it{#eta}_{MC}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    // h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]->GetXaxis()->SetRangeUser(-249,249);
    // h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]->GetYaxis()->SetRangeUser(-249,249);
    h_TMstudiesCalo_dEta_vsEta[icalo][_combCalo[icalo]]->GetXaxis()->SetRangeUser(-0.24,0.24);
    h_TMstudiesCalo_dEta_vsEta[icalo][_combCalo[icalo]]->Draw("col");
    //  TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
    //  el4->SetFillStyle(0);
    //  el4->SetLineColor(kGray+2);
    //  el4->SetLineWidth(2);
    //  el4->SetLineStyle(7);
    //  el4->Draw("same");
      // drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    float matchingwindoweta = ReturnEtaCaloMatching(icalo, _combCalo[icalo]);
    if(plotUncut){
      float minycutline = 0;
      float maxycutline = 0;
      if(GetCaloDirection(icalo)==0){
        minycutline = -4;
        maxycutline = -2;
      } else if(GetCaloDirection(icalo)==1){
        minycutline = -2;
        maxycutline = 1.4;
      } else if(GetCaloDirection(icalo)==2){
        minycutline = 1;
        maxycutline = 3.5;
      }
      DrawGammaLines(-matchingwindoweta, -matchingwindoweta, minycutline, maxycutline, 3, linecolx, 2);
      DrawGammaLines(matchingwindoweta, matchingwindoweta, minycutline, maxycutline,   3, linecolx, 2);
    }
      drawLatexAdd(Form("%s and %s #pi^{#pm} clusters matched",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),0.17,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cXYplot->Print(Form("%s/2D_caloMatching_deltaEta_vsEta_%s_%s.%s", outputDir.Data(),str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data(), suffix.Data()));
  }
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(!_curr_active_calo[icalo]) continue;
    if(!IsHCALCalorimeter(icalo)) continue;
    cout << "plotting calo " << str_calorimeterPlot[icalo].Data() << endl;
    bool isFwd = IsForwardCalorimeter(icalo);
    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudiesCalo_dPhi_vsEta[icalo][_combCalo[icalo]],Form("#it{#varphi}_{%s}-#it{#varphi}_{%s}",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()), "#it{#eta}_{MC}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    // h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]->GetXaxis()->SetRangeUser(-249,249);
    // h_TMstudiesCalo_2D_delta[icalo][_combCalo[icalo]]->GetYaxis()->SetRangeUser(-249,249);
    h_TMstudiesCalo_dPhi_vsEta[icalo][_combCalo[icalo]]->GetXaxis()->SetRangeUser(-0.24,0.24);
    h_TMstudiesCalo_dPhi_vsEta[icalo][_combCalo[icalo]]->Draw("col");
    //  TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
    //  el4->SetFillStyle(0);
    //  el4->SetLineColor(kGray+2);
    //  el4->SetLineWidth(2);
    //  el4->SetLineStyle(7);
    //  el4->Draw("same");
      // drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%s and %s #pi^{#pm} clusters matched",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),0.17,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    float matchingwindowphi = ReturnPhiCaloMatching(icalo, _combCalo[icalo]);
    if(plotUncut){
      float minycutline = 0;
      float maxycutline = 0;
      if(GetCaloDirection(icalo)==0){
        minycutline = -4;
        maxycutline = -2;
      } else if(GetCaloDirection(icalo)==1){
        minycutline = -2;
        maxycutline = 1.4;
      } else if(GetCaloDirection(icalo)==2){
        minycutline = 1;
        maxycutline = 3.5;
      }
      DrawGammaLines(-matchingwindowphi, -matchingwindowphi, minycutline, maxycutline, 3, linecolx, 2);
      DrawGammaLines(matchingwindowphi, matchingwindowphi, minycutline, maxycutline,   3, linecolx, 2);
    }

    cXYplot->Print(Form("%s/2D_caloMatching_deltaPhi_vsEta_%s_%s.%s", outputDir.Data(),str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data(), suffix.Data()));
  }



  TH1F*  h_TMstudies_clusSpec_Reb[_active_calo] = {nullptr};                        // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Reb_Matched[_active_calo] = {nullptr};                // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Reb_MatchedCalo[_active_calo] = {nullptr};            // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Reb_MatchedTrackAndCalo[_active_calo] = {nullptr};    // [calorimeter_enum][algorithm_enum]


  TH1F*  h_TMstudies_clusSpec_Cut_Reb[_active_calo] = {nullptr};                        // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Cut_Reb_Matched[_active_calo] = {nullptr};                // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Cut_Reb_MatchedCalo[_active_calo] = {nullptr};            // [calorimeter_enum]
  TH1F*  h_TMstudies_clusSpec_Cut_Reb_MatchedTrackAndCalo[_active_calo] = {nullptr};    // [calorimeter_enum][algorithm_enum]

  TH1F* h_CaloMatchingEffi_MatchedCalo_vsAll[_active_calo] = {nullptr};
  TH1F* h_CaloMatchingEffi_MatchedCalo_vsTM[_active_calo] = {nullptr};
  TH1F* h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[_active_calo] = {nullptr};
  TH1F* h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[_active_calo] = {nullptr};

  TH1F* h_CaloMatchingEffi_MatchedCalo_Cut_vsAll[_active_calo] = {nullptr};
  TH1F* h_CaloMatchingEffi_MatchedCalo_Cut_vsTM[_active_calo] = {nullptr};
  TH1F* h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll[_active_calo] = {nullptr};
  TH1F* h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM[_active_calo] = {nullptr};
  

  for(int icalo=0;icalo<_active_calo;icalo++){
    if(!_curr_active_calo[icalo]) continue;
    if(!IsHCALCalorimeter(icalo)) continue;
    h_TMstudies_clusSpec_Reb[icalo]                     = (TH1F*)h_TMstudies_clusSpec[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Reb_%s",str_calorimeterPlot[icalo].Data()), partP);
    h_TMstudies_clusSpec_Reb_Matched[icalo]             = (TH1F*)h_TMstudies_clusSpec_Matched[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Reb_Matched_%s",str_calorimeterPlot[icalo].Data()), partP);
    h_TMstudies_clusSpec_Reb_MatchedCalo[icalo]         = (TH1F*)h_TMstudies_clusSpec_MatchedCalo[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Reb_MatchedCalo_%s",str_calorimeterPlot[icalo].Data()), partP);
    h_TMstudies_clusSpec_Reb_MatchedTrackAndCalo[icalo] = (TH1F*)h_TMstudies_clusSpec_MatchedTrackAndCalo[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Reb_MatchedTrackAndCalo_%s",str_calorimeterPlot[icalo].Data()), partP);

    h_TMstudies_clusSpec_Cut_Reb[icalo]                     = (TH1F*)h_TMstudies_clusSpec_Cut[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Cut_Reb_%s",str_calorimeterPlot[icalo].Data()), partP);
    h_TMstudies_clusSpec_Cut_Reb_Matched[icalo]             = (TH1F*)h_TMstudies_clusSpec_Matched_Cut[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Cut_Reb_Matched_%s",str_calorimeterPlot[icalo].Data()), partP);
    h_TMstudies_clusSpec_Cut_Reb_MatchedCalo[icalo]         = (TH1F*)h_TMstudies_clusSpec_MatchedCalo_Cut[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Cut_Reb_MatchedCalo_%s",str_calorimeterPlot[icalo].Data()), partP);
    h_TMstudies_clusSpec_Cut_Reb_MatchedTrackAndCalo[icalo] = (TH1F*)h_TMstudies_clusSpec_MatchedTrackAndCalo_Cut[icalo]->Rebin(nP-2, Form("h_TMstudies_clusSpec_Cut_Reb_MatchedTrackAndCalo_%s",str_calorimeterPlot[icalo].Data()), partP);

    cout << "calculating efficiencies for " << str_calorimeterPlot[icalo].Data() << endl;
    h_CaloMatchingEffi_MatchedCalo_vsAll[icalo] = (TH1F*) h_TMstudies_clusSpec_Reb_MatchedCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedCalo_vsAll_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedCalo_vsAll[icalo]->Divide(h_CaloMatchingEffi_MatchedCalo_vsAll[icalo],h_TMstudies_clusSpec_Reb[icalo], 1, 1,"B" );
    h_CaloMatchingEffi_MatchedCalo_vsTM[icalo] = (TH1F*) h_TMstudies_clusSpec_Reb_MatchedCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedCalo_vsTM_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedCalo_vsTM[icalo]->Divide(h_CaloMatchingEffi_MatchedCalo_vsTM[icalo],h_TMstudies_clusSpec_Reb_Matched[icalo], 1, 1,"B" );

    h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo] = (TH1F*) h_TMstudies_clusSpec_Reb_MatchedTrackAndCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo]->Divide(h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo],h_TMstudies_clusSpec_Reb[icalo], 1, 1,"B" );
    h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo] = (TH1F*) h_TMstudies_clusSpec_Reb_MatchedTrackAndCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo]->Divide(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo],h_TMstudies_clusSpec_Reb_Matched[icalo], 1, 1,"B" );

    cout << "calculating efficiencies for " << str_calorimeterPlot[icalo].Data() << endl;
    h_CaloMatchingEffi_MatchedCalo_Cut_vsAll[icalo] = (TH1F*) h_TMstudies_clusSpec_Cut_Reb_MatchedCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedCalo_Cut_vsAll_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedCalo_Cut_vsAll[icalo]->Divide(h_CaloMatchingEffi_MatchedCalo_Cut_vsAll[icalo],h_TMstudies_clusSpec_Cut_Reb[icalo], 1, 1,"B" );
    h_CaloMatchingEffi_MatchedCalo_Cut_vsTM[icalo] = (TH1F*) h_TMstudies_clusSpec_Cut_Reb_MatchedCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedCalo_Cut_vsTM_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedCalo_Cut_vsTM[icalo]->Divide(h_CaloMatchingEffi_MatchedCalo_Cut_vsTM[icalo],h_TMstudies_clusSpec_Cut_Reb_Matched[icalo], 1, 1,"B" );

    h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll[icalo] = (TH1F*) h_TMstudies_clusSpec_Cut_Reb_MatchedTrackAndCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll[icalo]->Divide(h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll[icalo],h_TMstudies_clusSpec_Cut_Reb[icalo], 1, 1,"B" );
    h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM[icalo] = (TH1F*) h_TMstudies_clusSpec_Cut_Reb_MatchedTrackAndCalo[icalo]->Clone(Form("h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM_%s",str_calorimeterPlot[icalo].Data()));
    h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM[icalo]->Divide(h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM[icalo],h_TMstudies_clusSpec_Cut_Reb_Matched[icalo], 1, 1,"B" );
  }


  Style_t markerStyleCalo[maxcalo] = {20, 47, 29};
  Style_t markerStyleCalo3[maxcalo] = {24, 46, 30};
  Style_t markerStyleCalo2[maxcalo] = {21, 34, 43};
  Style_t markerStyleCalo4[maxcalo] = {25, 28, 42};
  Size_t markerSizeCalo[maxcalo] = {1.5, 1.8, 2.2};

  Color_t colorPID_light[nPID]          = {kGray+1, kRed-6, kGreen-6, kCyan-6, kBlue+1, kOrange};
  Color_t colorPID_light2[nPID]         = {kGray+2, kRed-2, kGreen-2, kCyan-6, kBlue+1, kOrange};
  Color_t colorPID_light3[nPID]         = {kGray+3, kRed-8, kGreen-8, kCyan-8, kBlue+1, kOrange};

  TCanvas* canvasMatchingEffi       = new TCanvas("canvasMatchingEffi", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasMatchingEffi,  0.1, 0.01, 0.015, 0.095);
  // canvasMatchingEffi->SetLogy(1);
  canvasMatchingEffi->SetLogx(1);
  TH2F * histoDummyMatchingEffi;
  histoDummyMatchingEffi                = new TH2F("histoDummyMatchingEffi", "histoDummyMatchingEffi",1000, 0.51, 39, 1000, 0.0, 1.24 );
  SetStyleHistoTH2ForGraphs( histoDummyMatchingEffi, "#it{E}_{clus,HCal} (GeV/#it{c})", "#it{#epsilon}_{calomatch}",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1, 510, 510);//(#times #epsilon_{pur})
  histoDummyMatchingEffi->GetYaxis()->SetLabelOffset(0.001);
  histoDummyMatchingEffi->GetXaxis()->SetNoExponent();
  histoDummyMatchingEffi->GetXaxis()->SetMoreLogLabels(kTRUE);

  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    if(!_curr_active_calo[icalo]) continue;
    if(!IsHCALCalorimeter(icalo)) continue;
  histoDummyMatchingEffi->DrawCopy();
  // legendPi0Merge1           = GetAndSetLegend2(0.50, 0.14, 0.70, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  // legendPi0Merge2           = GetAndSetLegend2(0.60, 0.14, 0.80, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  // legendPi0Merge3           = GetAndSetLegend2(0.70, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    TLegend* legendPi0Merge4           = GetAndSetLegend2(0.65, 0.14, 0.85, 0.14+(3.5*textSizeLabelsRel),1.2*textSizeLabelsPixel,1);
    if(icalo==kHCALIN)legendPi0Merge4           = GetAndSetLegend2(0.15, 0.54, 0.35, 0.54+(3.5*textSizeLabelsRel),1.2*textSizeLabelsPixel,1);
    int icaloPlot = 0;
    if(icalo==kLFHCAL) icaloPlot = 0;
    if(icalo==kHCALOUT) icaloPlot = 1;
    if(icalo==kHCALIN) icaloPlot = 2;
    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedCalo_vsAll[icalo], markerStyleCalo[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID[icaloPlot] , colorPID[icaloPlot]);
    h_CaloMatchingEffi_MatchedCalo_vsAll[icalo]->Draw("same");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_vsAll[icalo],"N_{clus,CM} / N_{clus}","p");

    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo], markerStyleCalo3[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light2[icaloPlot] , colorPID_light2[icaloPlot]);
    h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo]->Draw("same");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo],"N_{clus,CM,TM} / N_{clus}","p");

    // DrawGammaSetMarker(h_CaloMatchingEffi_MatchedCalo_vsTM[icalo], markerStyleCalo4[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light2[icaloPlot] , colorPID_light2[icaloPlot]);
    // h_CaloMatchingEffi_MatchedCalo_vsTM[icalo]->Draw("same");
    // legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_vsTM[icalo],"N_{clus,CM} / N_{clus,TM}","p");

    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo], markerStyleCalo2[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light[icaloPlot] , colorPID_light[icaloPlot]);
    h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo]->Draw("same");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo],"N_{clus,CM,TM} / N_{clus,TM}","p");

    legendPi0Merge4->Draw();
    float toplegval = 0.90;

    drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
    drawLatexAdd(Form("%s-%s cluster matching",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("1.6#cdot",0.55,0.35,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#sigma",0.52,0.30,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#scale[0.75]{FWHM}",0.58,0.30,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#epsilon_{95%}",0.70,0.30,textSizeLabelsRel,false,false,false);
    // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

    DrawGammaLines(0.51, 39, 1,1,2,kGray+1, 2);


    canvasMatchingEffi->Update();
    canvasMatchingEffi->Print(Form("%s/MatchingEffi_%s_%s.%s",outputDir.Data(),str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data(),suffix.Data()));

  histoDummyMatchingEffi->DrawCopy();
  // legendPi0Merge1           = GetAndSetLegend2(0.50, 0.14, 0.70, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  // legendPi0Merge2           = GetAndSetLegend2(0.60, 0.14, 0.80, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  // legendPi0Merge3           = GetAndSetLegend2(0.70, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    legendPi0Merge4           = GetAndSetLegend2(0.65, 0.14, 0.85, 0.14+(3.5*textSizeLabelsRel),1.2*textSizeLabelsPixel,1);
    if(icalo==kHCALIN)legendPi0Merge4           = GetAndSetLegend2(0.15, 0.54, 0.35, 0.54+(3.5*textSizeLabelsRel),1.2*textSizeLabelsPixel,1);
    icaloPlot = 0;
    if(icalo==kLFHCAL) icaloPlot = 0;
    if(icalo==kHCALOUT) icaloPlot = 1;
    if(icalo==kHCALIN) icaloPlot = 2;
    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedCalo_Cut_vsAll[icalo], markerStyleCalo[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID[icaloPlot] , colorPID[icaloPlot]);
    h_CaloMatchingEffi_MatchedCalo_Cut_vsAll[icalo]->Draw("same");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_Cut_vsAll[icalo],"N_{clus,CM} / N_{clus}","p");

    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll[icalo], markerStyleCalo3[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light2[icaloPlot] , colorPID_light2[icaloPlot]);
    h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll[icalo]->Draw("same");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsAll[icalo],"N_{clus,CM,TM} / N_{clus}","p");

    // DrawGammaSetMarker(h_CaloMatchingEffi_MatchedCalo_Cut_vsTM[icalo], markerStyleCalo4[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light2[icaloPlot] , colorPID_light2[icaloPlot]);
    // h_CaloMatchingEffi_MatchedCalo_Cut_vsTM[icalo]->Draw("same");
    // legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_Cut_vsTM[icalo],"N_{clus,CM} / N_{clus,TM}","p");

    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM[icalo], markerStyleCalo2[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light[icaloPlot] , colorPID_light[icaloPlot]);
    h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM[icalo]->Draw("same");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_Cut_vsTM[icalo],"N_{clus,CM,TM} / N_{clus,TM}","p");

    legendPi0Merge4->Draw();
    toplegval = 0.90;

    drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
    drawLatexAdd(Form("%s-%s cluster matching",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("1.6#cdot",0.55,0.35,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#sigma",0.52,0.30,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#scale[0.75]{FWHM}",0.58,0.30,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#epsilon_{95%}",0.70,0.30,textSizeLabelsRel,false,false,false);
    // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
    // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

    DrawGammaLines(0.51, 39, 1,1,2,kGray+1, 2);


    canvasMatchingEffi->Update();
    canvasMatchingEffi->Print(Form("%s/MatchingEffi_EtaCut_%s_%s.%s",outputDir.Data(),str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data(),suffix.Data()));
  }

  histoDummyMatchingEffi->DrawCopy();
  // legendPi0Merge1           = GetAndSetLegend2(0.50, 0.14, 0.70, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  // legendPi0Merge2           = GetAndSetLegend2(0.60, 0.14, 0.80, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
  // legendPi0Merge3           = GetAndSetLegend2(0.70, 0.14, 0.90, 0.14+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel,1);
    // TLegend* legendPi0Merge4           = GetAndSetLegend2(0.15, 0.74-(7.5*textSizeLabelsRel), 0.35, 0.74,1.1*textSizeLabelsPixel,1);
    TLegend* legendPi0Merge4           = GetAndSetLegend2(0.15, 0.74-(3.5*textSizeLabelsRel), 0.65, 0.74,1.1*textSizeLabelsPixel,2);
    legendPi0Merge4->SetMargin(0.1);
    // if(icalo==kHCALIN)legendPi0Merge4           = GetAndSetLegend2(0.15, 0.54, 0.35, 0.54+(3.5*textSizeLabelsRel),1.2*textSizeLabelsPixel,1);
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){

    if(!_curr_active_calo[icalo]) continue;
    if(!IsHCALCalorimeter(icalo)) continue;
    int icaloPlot = 0;
    if(icalo==kLFHCAL) icaloPlot = 0;
    if(icalo==kHCALOUT) icaloPlot = 1;
    if(icalo==kHCALIN) icaloPlot = 2;
    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedCalo_vsAll[icalo], markerStyleCalo[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID[icaloPlot] , colorPID[icaloPlot]);
    h_CaloMatchingEffi_MatchedCalo_vsAll[icalo]->Draw("same");
    // legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_vsAll[icalo],"N_{clus,CM} / N_{clus}","p");

    // DrawGammaSetMarker(h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo], markerStyleCalo3[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light2[icaloPlot] , colorPID_light2[icaloPlot]);
    // h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo]->Draw("same");
    // legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_vsAll[icalo],"N_{clus,CM,TM} / N_{clus}","p");

    // DrawGammaSetMarker(h_CaloMatchingEffi_MatchedCalo_vsTM[icalo], markerStyleCalo4[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light2[icaloPlot] , colorPID_light2[icaloPlot]);
    // h_CaloMatchingEffi_MatchedCalo_vsTM[icalo]->Draw("same");
    // legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_vsTM[icalo],"N_{clus,CM} / N_{clus,TM}","p");

    DrawGammaSetMarker(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo], markerStyleCalo3[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light2[icaloPlot] , colorPID_light2[icaloPlot]);
    // DrawGammaSetMarker(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo], markerStyleCalo2[icaloPlot], 1.5*markerSizeCalo[icaloPlot], colorPID_light[icaloPlot] , colorPID_light[icaloPlot]);
    h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo]->Draw("same");
    // legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo],"N_{clus,CM,TM} / N_{clus,TM}","p");
  }
  // legendPi0Merge4->AddEntry((TObject*)0,"N_{clus,CM} / N_{clus}"," ");
  // for (Int_t icalo = 0; icalo < maxcalo; icalo++){
  //   if(!_curr_active_calo[icalo]) continue;
  //   if(!IsHCALCalorimeter(icalo)) continue;
  //     legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_vsAll[icalo],Form("%s-%s",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),"p");
  // }
  // legendPi0Merge4->AddEntry((TObject*)0,"N_{clus,CM,TM} / N_{clus,TM}"," ");
  // for (Int_t icalo = 0; icalo < maxcalo; icalo++){
  //   if(!_curr_active_calo[icalo]) continue;
  //   if(!IsHCALCalorimeter(icalo)) continue;
  //   legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo],Form("%s-%s",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),"p");
  // }
  // legendPi0Merge4->AddEntry((TObject*)0," "," ");
  legendPi0Merge4->AddEntry((TObject*)0,"#it{#epsilon}_{CM,all}"," ");
  legendPi0Merge4->AddEntry((TObject*)0,"#it{#epsilon}_{CM,TM}"," ");
  // legendPi0Merge4->AddEntry((TObject*)0,"#frac{N_{clus,CM}}{N_{clus}}"," ");
  // legendPi0Merge4->AddEntry((TObject*)0,"#frac{N_{clus,CM,TM}}{N_{clus,TM}}"," ");
  // legendPi0Merge4->AddEntry((TObject*)0,"N_{clus,CM} / N_{clus}"," ");
  // legendPi0Merge4->AddEntry((TObject*)0,"N_{clus,CM,TM} / N_{clus,TM}"," ");
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    if(!_curr_active_calo[icalo]) continue;
    if(!IsHCALCalorimeter(icalo)) continue;
    // legendPi0Merge4->AddEntry((TObject*)0,Form("%s-%s",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data())," ");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedCalo_vsAll[icalo]," ","p");
    legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo]," ","p");
  }
  // legendPi0Merge4->AddEntry((TObject*)0,"N_{clus,CM,TM} / N_{clus,TM}"," ");
  // for (Int_t icalo = 0; icalo < maxcalo; icalo++){
  //   if(!_curr_active_calo[icalo]) continue;
  //   if(!IsHCALCalorimeter(icalo)) continue;
  //   legendPi0Merge4->AddEntry(h_CaloMatchingEffi_MatchedTrackAndCalo_vsTM[icalo],Form("%s-%s",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),"p");
  // }
  legendPi0Merge4->Draw();
  float toplegval = 0.90;

  drawLatexAdd(perfLabel.Data(),0.15,toplegval,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(Form("%s-%s cluster matching",str_calorimeterPlot[icalo].Data(),str_calorimeterPlot[_combCalo[icalo]].Data()),0.15,toplegval-0.05,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("1.6#cdot",0.55,0.35,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#sigma",0.52,0.30,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#scale[0.75]{FWHM}",0.58,0.30,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#epsilon_{95%}",0.70,0.30,textSizeLabelsRel,false,false,false);
  // drawLatexAdd(caloName[icalo].Data(),0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);
  // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,toplegval-0.1,textSizeLabelsRel,false,false,false);

  DrawGammaLines(0.51, 39, 1,1,2,kGray+1, 2);


  canvasMatchingEffi->Update();
  canvasMatchingEffi->Print(Form("%s/MatchingEffiCalos_Paper.%s",outputDir.Data(),suffix.Data()));

}

