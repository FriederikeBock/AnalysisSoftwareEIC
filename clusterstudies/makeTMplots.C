#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))


  const Int_t nEta                = 15;                                        
  Double_t partEta[nEta+1]        = { -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.2, -0.4, 0.4, 1.2, 
                                       1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

enum calotype {
    kFHCAL         = 0,
    kFEMC         = 1,
    kDRCALO        = 2,
    kEEMC         = 3,
    kCEMC         = 4,
    kEHCAL         = 5,
    kHCALIN       = 6,
    kHCALOUT       = 7,
    kLFHCAL        = 8,
    kEEMCG         = 9,
    kBECAL         = 10
};
const int maxcalo = 11;
int calogeomindex[maxcalo] = {0};

TString str_calorimeter[maxcalo] = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC", "EHCAL", "HCALIN", "HCALOUT", "LFHCAL", "EEMCG", "BECAL"};
const int _active_calo = 11;
bool _curr_active_calo[maxcalo] = {false, true, false, true, false, true, false, false, true, true, true};

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


bool IsForwardCalorimeter(int caloID){
  switch (caloID){
    case kDRCALO: return true;
    case kFHCAL: return true;
    case kFEMC: return true;
    case kEHCAL: return true;
    case kEEMC: return true;
    case kHCALIN: return false;
    case kHCALOUT: return false;
    case kCEMC: return false;
    case kEEMCG: return true;
    case kLFHCAL: return true;
    case kBECAL: return false;
    default:
      cout << "IsForwardCalorimeter: caloID " << caloID << " not defined, returning false" << endl;
      return false;
  }
  return false;
}


float ReturnTrackMatchingWindowForCalo(int caloID){
  switch (caloID){
    case kDRCALO: return 5;
    case kFHCAL: return 10;
    case kFEMC: return 5.5;
    case kEHCAL: return 10;
    case kEEMC: return 4.0;
    case kHCALIN: return 0.1;
    case kHCALOUT: return 0.1;
    case kCEMC: return 0.05;
    case kEEMCG: return 7;
    case kLFHCAL: return 10;
    case kBECAL: return 0.2;
    default:
      cout << "ReturnTrackMatchingWindowForCalo: caloID " << caloID << " not defined, returning -1" << endl;
      return -1;
  }
  return -1;
}
void makeTMplots(
    // TString nameFile  = "/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/NEWINPUT_TMSTUD_SINGLEP/output_TMSTUD.root",
    TString nameFile  = "/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/NEWINPUT_TMSTUD_SINGLEP_NEW/output_TMSTUD.root",
    TString suffix            = "pdf"
){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFile = new TFile(nameFile.Data());

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;

// if(h_TMstudies_x_y_track)h_TMstudies_x_y_track->Write();
//   for(int icalo=0;icalo<_active_calo;icalo++){
//     for(int icase=0;icase<diffCases_TMStudies;icase++){
//       if(h_TMstudies_isMatched[icase][icalo])h_TMstudies_isMatched[icase][icalo]->Write();
//       if(h_TMstudies_isMatched_E[icase][icalo])h_TMstudies_isMatched_E[icase][icalo]->Write();
//       if(h_TMstudies_isMatched_PDG[icase][icalo])h_TMstudies_isMatched_PDG[icase][icalo]->Write();
//       if(h_TMstudies_isMatched_findable[icase][icalo])h_TMstudies_isMatched_findable[icase][icalo]->Write();
//       if(h_TMstudies_isMatched_findable_E[icase][icalo])h_TMstudies_isMatched_findable_E[icase][icalo]->Write();
//       if(h_TMstudies_isMatched_findable_PDG[icase][icalo])h_TMstudies_isMatched_findable_PDG[icase][icalo]->Write();
//       if(h_TMstudies_x_y_projection[icalo]&&icase==0)h_TMstudies_x_y_projection[icalo]->Write();
//       for(int ialgo=0;ialgo<_active_algo;ialgo++){
//         if(h_TMstudies_x_y_clus[icalo][ialgo]&&icase==0) h_TMstudies_x_y_clus[icalo][ialgo]->Write();
//         if(h_TMstudies_2D_delta_track[icase][icalo][ialgo]) h_TMstudies_2D_delta_track[icase][icalo][ialgo]->Write();
//         if(h_TMstudies_2D_delta_projection[icase][icalo][ialgo]) h_TMstudies_2D_delta_projection[icase][icalo][ialgo]->Write();
//         if(h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo]) h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo]->Write();
//         if(h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo]) h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo]->Write();
//         if(h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo]) h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo]->Write();
//         if(h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo]) h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo]->Write();
//       }
//     }
//   }


  TH2F* h_TMstudies_x_y_clus[_active_calo] 	    = {NULL};
  TH2F* h_TMstudies_x_y_track[_active_calo] 	    = {NULL};
  TH2F* h_TMstudies_x_y_projection[_active_calo] 	    = {NULL};
  TH2F* h_TMstudies_all_2D_delta_track[_active_calo] 	    = {NULL};
  TH2F* h_TMstudies_all_2D_delta_projection[_active_calo] 	    = {NULL};

  TH1F* h_TMstudies_clusSpec_Eta[_active_calo][nEta+1] 	    = {{NULL}};
  TH1F* h_TMstudies_clusSpec_Matched_Eta[_active_calo][nEta+1] 	    = {{NULL}};
  TH1F* h_TMstudies_MatchingEffi_Eta[_active_calo][nEta+1] 	    = {{NULL}};

  TH1F* h_TMstudies_clusSpec[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_clusSpec_Matched[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_clusSpec_MatchedCalo[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_clusSpec_MatchedCalo_EDiff1[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_clusSpec_MatchedCalo_EDiff2[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_clusSpec_special[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_clusSpec_special_Matched[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_MatchingEffi[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_MatchingEffi_special[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_MatchingEffi_Calo[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_MatchingEffi_Calo_EDiff1[_active_calo] 	    = {NULL};
  TH1F* h_TMstudies_MatchingEffi_Calo_EDiff2[_active_calo] 	    = {NULL};

  for(int icalo=0;icalo<_active_calo;icalo++){
    bool isFwd = IsForwardCalorimeter(icalo);
    TString str2DhistoName = isFwd ? "x_y" : "eta_phi";

    h_TMstudies_x_y_clus[icalo] 	= (TH2F*) inputFile->Get(Form("%s/h_TMstudies_%s_clus_%s_MA",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()));
    if(!h_TMstudies_x_y_clus[icalo]) cout << Form("%s/h_TMstudies_%s_clus_%s_MA",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()) << " not found" << endl;

    h_TMstudies_x_y_track[icalo] 	= (TH2F*) inputFile->Get(Form("%s/h_TMstudies_%s_track_%s",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()));
    if(!h_TMstudies_x_y_track[icalo]) cout << Form("%s/h_TMstudies_%s_track_%s",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()) << " not found" << endl;

    h_TMstudies_x_y_projection[icalo] 	= (TH2F*) inputFile->Get(Form("%s/h_TMstudies_%s_projection_%s",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()));
    if(!h_TMstudies_x_y_projection[icalo]) cout << Form("%s/h_TMstudies_%s_projection_%s",str_calorimeter[icalo].Data(),str2DhistoName.Data(),str_calorimeter[icalo].Data()) << " not found" << endl;

    h_TMstudies_all_2D_delta_track[icalo] 	= (TH2F*) inputFile->Get(Form("%s/h_TMstudies_all_2D_delta_track_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    if(!h_TMstudies_all_2D_delta_track[icalo]) cout << Form("%s/h_TMstudies_all_2D_delta_track_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()) << " not found" << endl;

    h_TMstudies_all_2D_delta_projection[icalo] 	= (TH2F*) inputFile->Get(Form("%s/h_TMstudies_all_2D_delta_projection_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    if(!h_TMstudies_all_2D_delta_projection[icalo]) cout << Form("%s/h_TMstudies_all_2D_delta_projection_%s_MA",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()) << " not found" << endl;


    h_TMstudies_clusSpec[icalo] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_%s",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_clusSpec_Matched[icalo] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_Matched_%s",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_clusSpec_MatchedCalo[icalo] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_MatchedCalo_%s",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_MatchedCalo_EDiff1_%s",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_MatchedCalo_EDiff2_%s",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_clusSpec_special[icalo] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_special_%s",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));
    h_TMstudies_clusSpec_special_Matched[icalo] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_special_Matched_%s",str_calorimeter[icalo].Data(),str_calorimeter[icalo].Data()));

    h_TMstudies_clusSpec[icalo]->Rebin(2);
    h_TMstudies_clusSpec_special[icalo]->Rebin(2);
    h_TMstudies_clusSpec_Matched[icalo]->Rebin(2);
    h_TMstudies_clusSpec_MatchedCalo[icalo]->Rebin(2);
    h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]->Rebin(2);
    h_TMstudies_clusSpec_special_Matched[icalo]->Rebin(2);
    h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]->Rebin(2);

    h_TMstudies_MatchingEffi[icalo] 	= (TH1F*) h_TMstudies_clusSpec_Matched[icalo]->Clone(Form("h_TMstudies_MatchingEffi_%s",str_calorimeter[icalo].Data()));
    h_TMstudies_MatchingEffi[icalo]->Sumw2();
    h_TMstudies_MatchingEffi[icalo]->Divide(h_TMstudies_MatchingEffi[icalo],h_TMstudies_clusSpec[icalo],1,1,"b");

    h_TMstudies_MatchingEffi_Calo[icalo] 	= (TH1F*) h_TMstudies_clusSpec_MatchedCalo[icalo]->Clone(Form("h_TMstudies_MatchingEffi_Calo_%s",str_calorimeter[icalo].Data()));
    h_TMstudies_MatchingEffi_Calo[icalo]->Sumw2();
    h_TMstudies_MatchingEffi_Calo[icalo]->Divide(h_TMstudies_MatchingEffi_Calo[icalo],h_TMstudies_clusSpec[icalo],1,1,"b");
    h_TMstudies_MatchingEffi_Calo_EDiff1[icalo] 	= (TH1F*) h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]->Clone(Form("h_TMstudies_MatchingEffi_Calo_EDiff1_%s",str_calorimeter[icalo].Data()));
    h_TMstudies_MatchingEffi_Calo_EDiff1[icalo]->Sumw2();
    h_TMstudies_MatchingEffi_Calo_EDiff1[icalo]->Divide(h_TMstudies_MatchingEffi_Calo_EDiff1[icalo],h_TMstudies_clusSpec[icalo],1,1,"b");

    h_TMstudies_MatchingEffi_Calo_EDiff2[icalo] 	= (TH1F*) h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]->Clone(Form("h_TMstudies_MatchingEffi_Calo_EDiff2_%s",str_calorimeter[icalo].Data()));
    h_TMstudies_MatchingEffi_Calo_EDiff2[icalo]->Sumw2();
    h_TMstudies_MatchingEffi_Calo_EDiff2[icalo]->Divide(h_TMstudies_MatchingEffi_Calo_EDiff2[icalo],h_TMstudies_clusSpec[icalo],1,1,"b");

    h_TMstudies_MatchingEffi_special[icalo] 	= (TH1F*) h_TMstudies_clusSpec_special_Matched[icalo]->Clone(Form("h_TMstudies_MatchingEffi_special_%s",str_calorimeter[icalo].Data()));
    h_TMstudies_MatchingEffi_special[icalo]->Sumw2();
    h_TMstudies_MatchingEffi_special[icalo]->Divide(h_TMstudies_MatchingEffi_special[icalo],h_TMstudies_clusSpec_special[icalo],1,1,"b");
    for (Int_t et = 0; et<nEta+1; et++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (et < nEta){
        etaMin = partEta[et];
        etaMax = partEta[et+1];
      }
      h_TMstudies_clusSpec_Eta[icalo][et] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_Eta_%1.1f_%1.1f_%s",str_calorimeter[icalo].Data(),etaMin,etaMax,str_calorimeter[icalo].Data()));
      h_TMstudies_clusSpec_Eta[icalo][et]->Rebin(4);
      h_TMstudies_clusSpec_Matched_Eta[icalo][et] 	= (TH1F*) inputFile->Get(Form("%s/h_TMstudies_clusSpec_Matched_Eta_%1.1f_%1.1f_%s",str_calorimeter[icalo].Data(),etaMin,etaMax,str_calorimeter[icalo].Data()));
      h_TMstudies_clusSpec_Matched_Eta[icalo][et]->Rebin(4);
      h_TMstudies_MatchingEffi_Eta[icalo][et] 	= (TH1F*) h_TMstudies_clusSpec_Matched_Eta[icalo][et]->Clone(Form("h_TMstudies_MatchingEffi_Eta_%1.1f_%1.1f_%s",etaMin,etaMax,str_calorimeter[icalo].Data()));
      h_TMstudies_MatchingEffi_Eta[icalo][et]->Sumw2();
      h_TMstudies_MatchingEffi_Eta[icalo][et]->Divide(h_TMstudies_MatchingEffi_Eta[icalo][et],h_TMstudies_clusSpec_Eta[icalo][et],1,1,"b");
    }
  }

  // 2D PLOT
  TCanvas* cXYplot = new TCanvas("cXYplot","",0,0,1000,1000);
  DrawGammaCanvasSettings( cXYplot, 0.13, 0.015, 0.015, 0.1);
  // cXYplot->SetLogz();
  for(int icalo=0;icalo<_active_calo;icalo++){
    bool isFwd = IsForwardCalorimeter(icalo);
    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudies_x_y_track[icalo],isFwd ? "#it{x} [cm]" : "#eta", isFwd ? "#it{y} [cm]" : "#phi", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    // h_TMstudies_x_y_track[icalo]->GetXaxis()->SetRangeUser(-249,249);
    // h_TMstudies_x_y_track[icalo]->GetYaxis()->SetRangeUser(-249,249);
    h_TMstudies_x_y_track[icalo]->Draw("col");
    //  TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
    //  el4->SetFillStyle(0);
    //  el4->SetLineColor(kGray+2);
    //  el4->SetLineWidth(2);
    //  el4->SetLineStyle(7);
    //  el4->Draw("same");
      // drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#pi^{#pm} / e^{-} tracks prop. to %s surface",str_calorimeter[icalo].Data()),0.17,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cXYplot->Print(Form("%s/x_y_track_%s.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
    // cXYplot->Print(Form("%s/IEtaIPhiFHCAL_1D.%s", outputDir.Data(), suffix.Data()));
    SetStyleHistoTH2ForGraphs(h_TMstudies_x_y_projection[icalo],isFwd ? "#it{x} [cm]" : "#eta", isFwd ? "#it{y} [cm]" : "#phi", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    // h_TMstudies_x_y_track[icalo]->GetXaxis()->SetRangeUser(-249,249);
    // h_TMstudies_x_y_track[icalo]->GetYaxis()->SetRangeUser(-249,249);
    h_TMstudies_x_y_projection[icalo]->Draw("col");
    //  TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
    //  el4->SetFillStyle(0);
    //  el4->SetLineColor(kGray+2);
    //  el4->SetLineWidth(2);
    //  el4->SetLineStyle(7);
    //  el4->Draw("same");
      // drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#pi^{#pm} / e^{-} tracks projections on %s",str_calorimeter[icalo].Data()),0.17,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cXYplot->Print(Form("%s/x_y_projection_%s.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
    // cXYplot->Print(Form("%s/IEtaIPhiFHCAL_1D.%s", outputDir.Data(), suffix.Data()));

    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudies_x_y_clus[icalo],isFwd ? "#it{x} [cm]" : "#eta", isFwd ? "#it{y} [cm]" : "#phi", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    // h_TMstudies_x_y_track[icalo]->GetXaxis()->SetRangeUser(-249,249);
    // h_TMstudies_x_y_track[icalo]->GetYaxis()->SetRangeUser(-249,249);
    h_TMstudies_x_y_clus[icalo]->Draw("col");
    //  TEllipse *el4 = new TEllipse(granularity*27.0,granularity*27.0,granularity*27,granularity*27,0,360,0);
    //  el4->SetFillStyle(0);
    //  el4->SetLineColor(kGray+2);
    //  el4->SetLineWidth(2);
    //  el4->SetLineStyle(7);
    //  el4->Draw("same");
      // drawLatexAdd("FHCAL, e-p: 10#times 250 GeV",0.13,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("%s, #pi^{#pm} / e^{-} clusters",str_calorimeter[icalo].Data()),0.17,0.92,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#it{E}^{FHCAL}_{tot} = %.1f GeV", totEnergy),0.13,0.14,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cXYplot->Print(Form("%s/x_y_clusters_%s.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
    // cXYplot->Print(Form("%s/IEtaIPhiFHCAL_1D.%s", outputDir.Data(), suffix.Data()));

    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudies_all_2D_delta_track[icalo],isFwd ? "#it{x}^{trk}-#it{x}^{cls} [cm]" : "#eta^{trk}-#eta^{cls}", isFwd ? "#it{y}^{trk}-#it{y}^{cls} [cm]" : "#phi^{trk}-#phi^{cls}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    h_TMstudies_all_2D_delta_track[icalo]->Draw("col");
      drawLatexAdd(Form("%s, #pi^{#pm} / e^{-} clusters and tracks",str_calorimeter[icalo].Data()),0.05,0.02,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    float matchingwindow = ReturnTrackMatchingWindowForCalo(icalo);
    Color_t linecolx = kGray+1;
    DrawGammaLines(-matchingwindow, -matchingwindow, -matchingwindow, matchingwindow, 5, linecolx, 1);
    DrawGammaLines(matchingwindow, matchingwindow, -matchingwindow, matchingwindow,   5, linecolx, 1);
    DrawGammaLines(-matchingwindow, matchingwindow, -matchingwindow, -matchingwindow, 5, linecolx, 1);
    DrawGammaLines(-matchingwindow, matchingwindow, matchingwindow, matchingwindow,   5, linecolx, 1);
    cXYplot->Print(Form("%s/2D_delta_tracks_clusters_%s.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));

    // 1D PLOT
    SetStyleHistoTH2ForGraphs(h_TMstudies_all_2D_delta_projection[icalo],isFwd ? "#it{x}^{trk}-#it{x}^{cls} [cm]" : "#eta^{trk}-#eta^{cls}", isFwd ? "#it{y}^{trk}-#it{y}^{cls} [cm]" : "#phi^{trk}-#phi^{cls}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.3);
    h_TMstudies_all_2D_delta_projection[icalo]->Draw("col");
    drawLatexAdd(Form("%s, #pi^{#pm} / e^{-} clusters and projections",str_calorimeter[icalo].Data()),0.05,0.02,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    DrawGammaLines(-matchingwindow, -matchingwindow, -matchingwindow, matchingwindow, 5, linecolx, 1);
    DrawGammaLines(matchingwindow, matchingwindow, -matchingwindow, matchingwindow,   5, linecolx, 1);
    DrawGammaLines(-matchingwindow, matchingwindow, -matchingwindow, -matchingwindow, 5, linecolx, 1);
    DrawGammaLines(-matchingwindow, matchingwindow, matchingwindow, matchingwindow,   5, linecolx, 1);
    cXYplot->Print(Form("%s/2D_delta_projection_clusters_%s.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
  }


  TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
  // canvasWeights->SetLogx();

  TH2F * histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0,50,1000,-0.01,1.1);
  SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{E}_{clus} (GeV)","#epsilon_{ TM}",0.04,0.045, 0.04,0.045, 0.9,0.8);
  histo2DWeights->GetXaxis()->SetMoreLogLabels();
  histo2DWeights->GetXaxis()->SetNoExponent();
  canvasWeights->cd();
  histo2DWeights->Draw("copy");
  // Color_t colorCalo[2] = {kRed+2,kGreen+2};
  // Style_t markerCalo[2] = {24,27};
    Color_t colorCalo[11]                          = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2, kOrange-2};
    Marker_t markerCalo    [11]                = {20, 20, 21, 34, 29, 33, 21, 27, 28, 30,46 };

  TLegend* legendInvY5TtionPaper    = GetAndSetLegend2(0.17, 0.18, 0.5, 0.18+0.05*2, 1.2*textSizeLabelsPixel);
  legendInvY5TtionPaper->SetMargin(0.2);
  if(_active_calo==11){
    legendInvY5TtionPaper    = GetAndSetLegend2(0.17, 0.18, 0.95, 0.18+0.05*2, 1.2*textSizeLabelsPixel);
    legendInvY5TtionPaper->SetMargin(0.2);
    legendInvY5TtionPaper->SetNColumns(3);
  }
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(!_curr_active_calo[icalo])continue;
    DrawGammaSetMarker(h_TMstudies_MatchingEffi[icalo], markerCalo[icalo], 2, colorCalo[icalo] , colorCalo[icalo]);
    h_TMstudies_MatchingEffi[icalo]->Draw("same");
    legendInvY5TtionPaper->AddEntry(h_TMstudies_MatchingEffi[icalo],Form("%s clusters",str_calorimeter[icalo].Data()),"p");
  }
  DrawGammaLines(0, 50, 1, 1,   3, kGray+2, 7);

  legendInvY5TtionPaper->Draw();

  canvasWeights->Print(Form("%s/trackmatchingeffi.%s", outputDir.Data(), suffix.Data()));

  histo2DWeights->Draw("copy");

  legendInvY5TtionPaper    = GetAndSetLegend2(0.17, 0.18, 0.5, 0.18+0.05*2, 1.2*textSizeLabelsPixel);
  legendInvY5TtionPaper->SetMargin(0.2);
  for(int icalo=0;icalo<_active_calo;icalo++){
    DrawGammaSetMarker(h_TMstudies_MatchingEffi_special[icalo], markerCalo[icalo], 2, colorCalo[icalo] , colorCalo[icalo]);
    h_TMstudies_MatchingEffi_special[icalo]->Draw("same");
    legendInvY5TtionPaper->AddEntry(h_TMstudies_MatchingEffi_special[icalo],Form("%s clusters",str_calorimeter[icalo].Data()),"p");
  }
  DrawGammaLines(0, 50, 1, 1,   3, kGray+2, 7);

  legendInvY5TtionPaper->Draw();

  canvasWeights->Print(Form("%s/trackmatchingeffi_special.%s", outputDir.Data(), suffix.Data()));

  histo2DWeights->Draw("copy");

  legendInvY5TtionPaper    = GetAndSetLegend2(0.17, 0.58, 0.5, 0.58+0.05*2, 1.2*textSizeLabelsPixel);
  legendInvY5TtionPaper->SetMargin(0.2);
  for(int icalo=0;icalo<_active_calo;icalo++){
    DrawGammaSetMarker(h_TMstudies_MatchingEffi_special[icalo], markerCalo[icalo], 2, colorCalo[icalo] , colorCalo[icalo]);
    h_TMstudies_MatchingEffi_special[icalo]->Draw("same");
    legendInvY5TtionPaper->AddEntry(h_TMstudies_MatchingEffi_special[icalo],Form("%s clusters",str_calorimeter[icalo].Data()),"p");
  }
  DrawGammaSetMarker(h_TMstudies_MatchingEffi_Calo[kFHCAL], 28, 2, kBlue+2 , kBlue+2);
  h_TMstudies_MatchingEffi_Calo[kFHCAL]->Draw("same");
  DrawGammaSetMarker(h_TMstudies_MatchingEffi_Calo_EDiff1[kFHCAL], 28, 2, kBlue-8 , kBlue-8);
  h_TMstudies_MatchingEffi_Calo_EDiff1[kFHCAL]->Draw("same");
  DrawGammaSetMarker(h_TMstudies_MatchingEffi_Calo_EDiff2[kFHCAL], 28, 2, kBlue+4 , kBlue+4); // FHCAL cluster > FEMC cluster
  h_TMstudies_MatchingEffi_Calo_EDiff2[kFHCAL]->Draw("same");
  DrawGammaLines(0, 50, 1, 1,   3, kGray+2, 7);

  legendInvY5TtionPaper->Draw();

  canvasWeights->Print(Form("%s/trackmatchingeffi_special2.%s", outputDir.Data(), suffix.Data()));
    
  for(int icalo=0;icalo<_active_calo;icalo++){

    histo2DWeights->Draw("copy");
    Color_t colorTrigg[10]                          = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};
    Marker_t markerTrigg    [10]                = {20, 20, 21, 34, 29, 33, 21, 27, 28, 30 };
    if(icalo==kFHCAL)
    legendInvY5TtionPaper    = GetAndSetLegend2(0.17, 0.50, 0.5, 0.50+0.05*nEta/2, 1.2*textSizeLabelsPixel);
    else
    legendInvY5TtionPaper    = GetAndSetLegend2(0.17, 0.18, 0.5, 0.18+0.05*nEta/2, 1.2*textSizeLabelsPixel);
    legendInvY5TtionPaper->SetMargin(0.2);
    // legendInvY5TtionPaper->SetNColumns(2);
    int currplot = 0;
    for (Int_t et = 0; et<nEta+1; et++){
      Double_t etaMin = partEta[0];
      Double_t etaMax = partEta[nEta];
      if (et < nEta){
        etaMin = partEta[et];
        etaMax = partEta[et+1];
      }
      if((icalo==kLFHCAL || icalo==kFEMC || icalo==kDRCALO || icalo==kFHCAL) && etaMin<1.2) continue;
      if((icalo==kEEMC || icalo==kEEMCG || icalo==kEHCAL ) && etaMin>=-1.2) continue;
      if((icalo==kCEMC || icalo==kBECAL || icalo==kHCALIN || icalo==kHCALOUT ) && (etaMin<-1.2 || etaMax>1.2)) continue;
      // DrawGammaSetMarker(h_TMstudies_MatchingEffi_Eta[icalo][et], markerCalo[icalo], 2, colorCalo[icalo] , colorCalo[icalo]);
      DrawGammaSetMarker(h_TMstudies_MatchingEffi_Eta[icalo][et], markerTrigg[currplot], 2, colorTrigg[currplot] , colorTrigg[currplot]);
      h_TMstudies_MatchingEffi_Eta[icalo][et]->Draw("same");
      legendInvY5TtionPaper->AddEntry(h_TMstudies_MatchingEffi_Eta[icalo][et],Form("%1.1f < #eta < %1.1f",etaMin,etaMax),"p");
      currplot++;
    }
      drawLatexAdd(Form("%s",str_calorimeter[icalo].Data()),0.85,0.13,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    legendInvY5TtionPaper->Draw();
  DrawGammaLines(0, 50, 1, 1,   3, kGray+2, 7);

    canvasWeights->Print(Form("%s/trackmatchingeffi_Etadep_%s.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
  }
}

