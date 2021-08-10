// #include <Riostream.h>
// #include "TMath.h"
// #include <stdlib.h>
// #include <fstream>
// #include <math.h>
// #include <TROOT.h>
// #include <TApplication.h>
// #include <TPaveLabel.h>
// #include <TSystem.h>
// #include <TFrame.h>
// #include <TStyle.h>
// #include <TString.h>
// #include "TGaxis.h"
// #include "TFile.h"
// #include "TH1F.h"
// #include "TH1D.h"
// #include "TH2F.h"
// #include "TF1.h"
// #include "TVirtualFitter.h"
// #include "TObject.h"
// #include "TCanvas.h"
// #include "TMultiGraph.h"
// #include "TLegend.h"
// #include "TDatabasePDG.h"
// #include "TMinuit.h"
// #include "TBenchmark.h"
// #include "TRandom.h"
// #include "TLatex.h"
// #include "TASImage.h"
// #include "TPostScript.h"
// #include "TGraphErrors.h"
// #include "TArrow.h"
// #include "TGraphAsymmErrors.h"
// #include "TGaxis.h"
// #include "TMarker.h"
// #include "Math/WrappedTF1.h"
// #include "Math/BrentRootFinder.h"
// #include "TFitResultPtr.h"
// #include "TFitResult.h"
#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value,const Float_t min_towers, const bool normalize, Float_t etalow, Float_t etahigh, Float_t calibconst = 1);
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);


void resolutionHCAL(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plotsResolutionHCAL/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Double_t energySingleCalibPlot = 8;
  //************************** Read data **************************************************
  const Int_t nEne = 25;
  const static Double_t partE[]   = {1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0,  10., 11., 13., 15., 20., 25., 30., 40., 50., 60., 70., 80., 150.};
  Int_t useE[nEne]                = {  0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1};
  Int_t useEgamma[nEne]           = {  0,   0,   1,   1,   1,   1,   0,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,   1};
  Double_t fitrange_low[nEne]     = {0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6,0.65, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7};
  Double_t fitrange_hig[nEne]     = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0};
  TH1F* h_pion_fhcal_E[nEne] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta1[nEne] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta2[nEne] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta3[nEne] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta4[nEne] 	    = {NULL};
  TH1F* h_pion_fhcal_E_allclus[nEne] 	    = {NULL};
  TH1F* h_gamma_fhcal_E[nEne] 	    = {NULL};
  Int_t activeE = 0;
  Int_t activeEgamma = 0;

  Int_t useEBothCPion[nEne]       = {  0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1};
  Int_t useEBothgamma[nEne]       = {  0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   1,   0};
  Int_t useEBothproton[nEne]      = {  0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0};
  TH1F* h_pion__dEta_hcal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dPhi_hcal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dR_hcal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dEta_ecal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dPhi_ecal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dR_ecal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dEta_ecalfhcal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dPhi_ecalfhcal_E[nEne] 	    = {NULL};
  TH1F* h_pion__dR_ecalfhcal_E[nEne] 	    = {NULL};
  TH1F* h_pion_ecalfhcal_E[nEne] 	    = {NULL};
  TH1F* h_pion_ecal_E[nEne] 	    = {NULL};
  TH1F* h_pion_ecalfhcalcomb_E[nEne] 	    = {NULL};
  TH1F* h_pion_ecalfhcalcomb_E_calib[nEne] 	    = {NULL};
  TH1F* h_pion_ecalfhcalcomb_E_calibrated[nEne] 	    = {NULL};
  TH2F* h2D_pion_ecalfhcal_E[nEne] 	    = {NULL};
  TH2F* h2D_NTowers_HCAL_pi[nEne] 	    = {NULL};
  TH2F* h2D_NTowers_HCAL_gamma[nEne] 	    = {NULL};
  TH2F* h2D_NTowers_HCAL_proton[nEne] 	    = {NULL};
  TH2F* h2D_NTowers_ECAL_pi[nEne] 	    = {NULL};
  TH2F* h2D_NTowers_ECAL_gamma[nEne] 	    = {NULL};
  TH2F* h2D_NTowers_ECAL_proton[nEne] 	    = {NULL};
  TH1F* h_gamma_ecalfhcal_E[nEne] 	    = {NULL};
  TH1F* h_gamma_ecal_E[nEne] 	    = {NULL};
  TH1F* h_gamma_ecal_E_calibrated[nEne] 	    = {NULL};
  TH1F* h_gamma_ecalfhcalcomb_E[nEne] 	    = {NULL};
  TH1F* h_proton_ecalfhcal_E[nEne] 	    = {NULL};
  TH1F* h_proton_ecal_E[nEne] 	    = {NULL};
  TH1F* h_proton_ecalfhcalcomb_E[nEne] 	    = {NULL};
  TH1F* h_proton_ecalfhcalcomb_E_calibrated[nEne] 	    = {NULL};
  TH2F* h2D_gamma_ecalfhcal_E[nEne] 	    = {NULL};
  TH2F* h2D_proton_ecalfhcal_E[nEne] 	    = {NULL};

  Int_t activeEboth = 0;
  Int_t activeEboth_gamma = 0;
  Int_t activeEboth_proton = 0;

  Double_t etalowdef = 1.3;
  Double_t etahighdef = 3.0;

  Double_t calibconstECAL = 8.36627e-01;
  Double_t paramsHCALcalib[5]= {2.65428e+00,4.17288e-01,1.91848e+00,-5.35529e+02,6.82686e+02}; //shift
  // Double_t paramsHCALcalib[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
  TF1* FunctionHCAL_calib   = new TF1("FunctionHCAL_calib",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   2, 220);
  FunctionHCAL_calib->SetParameters(paramsHCALcalib);

  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    cout << partE[epart] << endl;
    if(useE[epart]){
      TTree *const t_pion_fhcal =  (TTree *) (new TFile(Form("/home/nschmidt/fun4all_ForwardAna/EICdetector/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower"); //ntp_gshower
      if(t_pion_fhcal){
        activeE++;
        h_pion_fhcal_E[epart] 	= new TH1F(Form("h_pion_fhcal_E_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta1[epart] 	= new TH1F(Form("h_pion_fhcal_E_eta1_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta2[epart] 	= new TH1F(Form("h_pion_fhcal_E_eta2_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta3[epart] 	= new TH1F(Form("h_pion_fhcal_E_eta3_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta4[epart] 	= new TH1F(Form("h_pion_fhcal_E_eta4_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_allclus[epart] 	= new TH1F(Form("h_pion_fhcal_E_allclus_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_fhcal_E[epart], t_pion_fhcal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
        fill_histogram(h_pion_fhcal_E_eta1[epart], t_pion_fhcal, partE[epart], 0.3, 2, true, 1, 1.7);
        fill_histogram(h_pion_fhcal_E_eta2[epart], t_pion_fhcal, partE[epart], 0.3, 2, true, 1.7, 2.5);
        fill_histogram(h_pion_fhcal_E_eta3[epart], t_pion_fhcal, partE[epart], 0.3, 2, true, 2.5, 3.0);
        fill_histogram(h_pion_fhcal_E_eta4[epart], t_pion_fhcal, partE[epart], 0.3, 2, true, 3.0, 4.0);
        fill_histogram(h_pion_fhcal_E_allclus[epart], t_pion_fhcal, partE[epart], 0.3, 1, true, etalowdef, etahighdef);
      }
    }
    if(useEgamma[epart]){

      TTree *const t_gamma_fhcal =  (TTree *) (new TFile(Form("/home/nschmidt/fun4all_ForwardAna/FHCAL_gamma/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      if(t_gamma_fhcal){
        activeEgamma++;
        h_gamma_fhcal_E[epart] 	= new TH1F(Form("h_gamma_fhcal_E_%d",epart), "", 200, 0.0, 2.0);

        fill_histogram(h_gamma_fhcal_E[epart], t_gamma_fhcal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
    }
    if(useEBothCPion[epart]){
      TTree *const t_pion_fecalfhcal =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_cpion/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_pion_fecal =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_cpion/%1.1f/G4EICDetector_g4femc_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_pion_fecalfhcal_hcal =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_cpion/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_pion_fecalfhcal_ecal =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_cpion/%1.1f/G4EICDetector_g4femc_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      if(t_pion_fecalfhcal){
        activeEboth++;
        h_pion_ecalfhcal_E[epart] 	= new TH1F(Form("h_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_ecalfhcal_E[epart], t_pion_fecalfhcal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_pion_fecal){
        h_pion_ecal_E[epart] 	= new TH1F(Form("h_pion_ecal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_ecal_E[epart], t_pion_fecal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_pion_fecalfhcal_hcal && t_pion_fecalfhcal_ecal){
        h_pion__dEta_ecalfhcal_E[epart] 	= new TH1F(Form("h_pion__dEta_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_ecalfhcal_E[epart] 	= new TH1F(Form("h_pion__dPhi_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_ecalfhcal_E[epart] 	= new TH1F(Form("h_pion__dR_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dEta_hcal_E[epart] 	= new TH1F(Form("h_pion__dEta_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_hcal_E[epart] 	= new TH1F(Form("h_pion__dPhi_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_hcal_E[epart] 	= new TH1F(Form("h_pion__dR_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dEta_ecal_E[epart] 	= new TH1F(Form("h_pion__dEta_ecal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_ecal_E[epart] 	= new TH1F(Form("h_pion__dPhi_ecal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_ecal_E[epart] 	= new TH1F(Form("h_pion__dR_ecal_E_%d",epart), "", 200, -0.5, 0.5);


        h2D_NTowers_HCAL_pi[epart] 	= new TH2F(Form("h2D_NTowers_HCAL_pi%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_NTowers_ECAL_pi[epart] 	= new TH2F(Form("h2D_NTowers_ECAL_pi%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_pion_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0,partE[epart]*1.15, 100, 0,partE[epart]*1.15);
        h_pion_ecalfhcalcomb_E[epart] 	= new TH1F(Form("h_pion_ecalfhcalcomb_E_%d",epart), "", 100, 0.0, 2.0);
        h_pion_ecalfhcalcomb_E_calib[epart] 	= new TH1F(Form("h_pion_ecalfhcalcomb_E_calib_%d",epart), "", 100, 0.0, 2.0);
        h_pion_ecalfhcalcomb_E_calibrated[epart] 	= new TH1F(Form("h_pion_ecalfhcalcomb_E_calibrated_%d",epart), "", 100, 0.0, 2.0);
        Float_t measured_energy1,measured_energy2;
        t_pion_fecalfhcal_hcal->SetBranchAddress("e", &measured_energy1);
        t_pion_fecalfhcal_ecal->SetBranchAddress("e", &measured_energy2);
        Float_t true_eta1,true_eta2;
        t_pion_fecalfhcal_hcal->SetBranchAddress("geta", &true_eta1);
        t_pion_fecalfhcal_ecal->SetBranchAddress("geta", &true_eta2);
        Float_t true_phi1,true_phi2;
        t_pion_fecalfhcal_hcal->SetBranchAddress("gphi", &true_phi1);
        t_pion_fecalfhcal_ecal->SetBranchAddress("gphi", &true_phi2);
        Float_t rec_eta1,rec_eta2;
        t_pion_fecalfhcal_hcal->SetBranchAddress("eta", &rec_eta1);
        t_pion_fecalfhcal_ecal->SetBranchAddress("eta", &rec_eta2);
        Float_t rec_phi1,rec_phi2;
        t_pion_fecalfhcal_hcal->SetBranchAddress("phi", &rec_phi1);
        t_pion_fecalfhcal_ecal->SetBranchAddress("phi", &rec_phi2);
        Float_t nTowers1,nTowers2;
        t_pion_fecalfhcal_hcal->SetBranchAddress("ntowers", &nTowers1);
        t_pion_fecalfhcal_ecal->SetBranchAddress("ntowers", &nTowers2);

        Int_t nentries = Int_t(t_pion_fecalfhcal_hcal->GetEntries()) > Int_t(t_pion_fecalfhcal_ecal->GetEntries()) ? Int_t(t_pion_fecalfhcal_ecal->GetEntries()) : Int_t(t_pion_fecalfhcal_hcal->GetEntries());
        Float_t lastEvt = 0;
        for (Int_t i = 0; i < nentries; ++i) {
          if (t_pion_fecalfhcal_hcal->LoadTree(i) < 0 || t_pion_fecalfhcal_ecal->LoadTree(i) < 0)
            break;

          t_pion_fecalfhcal_hcal->GetEntry(i);
          t_pion_fecalfhcal_ecal->GetEntry(i);
          if (((true_eta1 > 1.3 && true_eta1 < 3.0)
            || (true_eta2 > 1.3 && true_eta2 < 3.0))  //1.7 < 3
          && (measured_energy1 > 0.35) // above MIP
          && (measured_energy1 > 0.35) // above MIP
          && nTowers1>=2 && nTowers2>=2 // actual cluster with >1 tower
          ){
            h2D_NTowers_HCAL_pi[epart]->Fill(measured_energy1, nTowers1);
            h2D_NTowers_ECAL_pi[epart]->Fill(measured_energy2, nTowers2);
            h2D_pion_ecalfhcal_E[epart]->Fill(measured_energy2, measured_energy1);
            h_pion_ecalfhcalcomb_E[epart]->Fill((measured_energy2+measured_energy1)/partE[epart]);
            h_pion_ecalfhcalcomb_E_calib[epart]->Fill(((measured_energy2/calibconstECAL)+measured_energy1)/partE[epart]);
            h_pion_ecalfhcalcomb_E_calibrated[epart]->Fill(((measured_energy2/calibconstECAL)+(measured_energy1/FunctionHCAL_calib->Eval(measured_energy1)))/partE[epart]);
            h_pion__dEta_ecalfhcal_E[epart]->Fill(rec_eta1-rec_eta2);
            h_pion__dPhi_ecalfhcal_E[epart]->Fill(rec_phi1-rec_phi2);
            h_pion__dR_ecalfhcal_E[epart]->Fill(TMath::Sqrt(TMath::Power(rec_phi1-rec_phi2,2)+TMath::Power(rec_eta1-rec_eta2,2)));
            h_pion__dEta_hcal_E[epart]->Fill(rec_eta1-true_eta1);
            h_pion__dPhi_hcal_E[epart]->Fill(rec_phi1-true_phi1);
            h_pion__dR_hcal_E[epart]->Fill(TMath::Sqrt(TMath::Power(rec_phi1-true_phi1,2)+TMath::Power(rec_eta1-true_eta1,2)));
            h_pion__dEta_ecal_E[epart]->Fill(rec_eta2-true_eta2);
            h_pion__dPhi_ecal_E[epart]->Fill(rec_phi2-true_phi2);
            h_pion__dR_ecal_E[epart]->Fill(TMath::Sqrt(TMath::Power(rec_phi2-true_phi2,2)+TMath::Power(rec_eta2-true_eta2,2)));
          }
        }
        h_pion_ecalfhcalcomb_E[epart]->Scale(1 / h_pion_ecalfhcalcomb_E[epart]->GetEntries());
        h_pion_ecalfhcalcomb_E_calib[epart]->Scale(1 / h_pion_ecalfhcalcomb_E_calib[epart]->GetEntries());
        h_pion_ecalfhcalcomb_E_calibrated[epart]->Scale(1 / h_pion_ecalfhcalcomb_E_calibrated[epart]->GetEntries());
        h_pion__dEta_ecalfhcal_E[epart]->Scale(1 / h_pion__dEta_ecalfhcal_E[epart]->GetEntries());
        h_pion__dPhi_ecalfhcal_E[epart]->Scale(1 / h_pion__dPhi_ecalfhcal_E[epart]->GetEntries());
        h_pion__dR_ecalfhcal_E[epart]->Scale(1 / h_pion__dR_ecalfhcal_E[epart]->GetEntries());
        h_pion__dEta_hcal_E[epart]->Scale(1 / h_pion__dEta_hcal_E[epart]->GetEntries());
        h_pion__dPhi_hcal_E[epart]->Scale(1 / h_pion__dPhi_hcal_E[epart]->GetEntries());
        h_pion__dR_hcal_E[epart]->Scale(1 / h_pion__dR_hcal_E[epart]->GetEntries());
        h_pion__dEta_ecal_E[epart]->Scale(1 / h_pion__dEta_ecal_E[epart]->GetEntries());
        h_pion__dPhi_ecal_E[epart]->Scale(1 / h_pion__dPhi_ecal_E[epart]->GetEntries());
        h_pion__dR_ecal_E[epart]->Scale(1 / h_pion__dR_ecal_E[epart]->GetEntries());
      }
    }
    if(useEBothgamma[epart]){
      TTree *const t_gamma_fecalfhcal_hcal1 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_gamma/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_gamma_fecalfhcal_ecal1 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_gamma/%1.1f/G4EICDetector_g4femc_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_gamma_fecalfhcal_hcal2 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_gamma/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_gamma_fecalfhcal_ecal2 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_gamma/%1.1f/G4EICDetector_g4femc_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      if(t_gamma_fecalfhcal_hcal1){
        activeEboth_gamma++;
        h_gamma_ecalfhcal_E[epart] 	= new TH1F(Form("h_gamma_ecalfhcal_E_%d",epart), "", 200, 0.0, 2.0);

        fill_histogram(h_gamma_ecalfhcal_E[epart], t_gamma_fecalfhcal_hcal1, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_gamma_fecalfhcal_ecal1){
        h_gamma_ecal_E[epart] 	= new TH1F(Form("h_gamma_ecal_E_%d",epart), "", 200, 0.0, 2.0);
        h_gamma_ecal_E_calibrated[epart] 	= new TH1F(Form("h_gamma_ecal_E_calibrated_%d",epart), "", 200, 0.0, 2.0);

        fill_histogram(h_gamma_ecal_E[epart], t_gamma_fecalfhcal_ecal1, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
        fill_histogram(h_gamma_ecal_E_calibrated[epart], t_gamma_fecalfhcal_ecal1, partE[epart], 0.3, 2, true, etalowdef, etahighdef,calibconstECAL);
      }
      if(t_gamma_fecalfhcal_hcal2 && t_gamma_fecalfhcal_ecal2){
        h2D_NTowers_HCAL_gamma[epart] 	= new TH2F(Form("h2D_NTowers_HCAL_gamma%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_NTowers_ECAL_gamma[epart] 	= new TH2F(Form("h2D_NTowers_ECAL_gamma%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_gamma_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_gamma_ecalfhcal_E_%d",epart), "", 100, 0,partE[epart]*1.15, 100, 0,partE[epart]*1.15);
        h_gamma_ecalfhcalcomb_E[epart] 	= new TH1F(Form("h_gamma_ecalfhcalcomb_E_%d",epart), "", 600, 0.0, 2.0);
        Float_t measured_energy1,measured_energy2;
        t_gamma_fecalfhcal_hcal2->SetBranchAddress("e", &measured_energy1);
        t_gamma_fecalfhcal_ecal2->SetBranchAddress("e", &measured_energy2);
        Float_t true_eta1,true_eta2;
        t_gamma_fecalfhcal_hcal2->SetBranchAddress("geta", &true_eta1);
        t_gamma_fecalfhcal_ecal2->SetBranchAddress("geta", &true_eta2);
        Float_t nTowers1,nTowers2;
        t_gamma_fecalfhcal_hcal2->SetBranchAddress("ntowers", &nTowers1);
        t_gamma_fecalfhcal_ecal2->SetBranchAddress("ntowers", &nTowers2);

        Int_t nentries = Int_t(t_gamma_fecalfhcal_hcal2->GetEntries()) > Int_t(t_gamma_fecalfhcal_ecal2->GetEntries()) ? Int_t(t_gamma_fecalfhcal_ecal2->GetEntries()) : Int_t(t_gamma_fecalfhcal_hcal2->GetEntries());
        Float_t lastEvt = 0;
        for (Int_t i = 0; i < nentries; ++i) {
          if (t_gamma_fecalfhcal_hcal2->LoadTree(i) < 0 || t_gamma_fecalfhcal_ecal2->LoadTree(i) < 0)
            break;

          t_gamma_fecalfhcal_hcal2->GetEntry(i);
          t_gamma_fecalfhcal_ecal2->GetEntry(i);
          if (((true_eta1 > 1.7 && true_eta1 < 3.0)
            || (true_eta2 > 1.7 && true_eta2 < 3.0))  //1.7 < 3
          && (measured_energy1 > 0.1) // above MIP
          && (measured_energy1 > 0.1) // above MIP
          // && nTowers1>=2 && nTowers2>=2 // actual cluster with >1 tower
          ){
            h2D_NTowers_HCAL_gamma[epart]->Fill(measured_energy1, nTowers1);
            h2D_NTowers_HCAL_gamma[epart]->Fill(measured_energy2, nTowers2);
            h2D_gamma_ecalfhcal_E[epart]->Fill(measured_energy2, measured_energy1);
            h_gamma_ecalfhcalcomb_E[epart]->Fill((measured_energy2+measured_energy1)/partE[epart]);
          }
        }
        h_gamma_ecalfhcalcomb_E[epart]->Scale(1 / h_gamma_ecalfhcalcomb_E[epart]->GetEntries());
      }
    }
    if(useEBothproton[epart]){
      TTree *const t_proton_fecalfhcal_hcal1 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_proton/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_proton_fecalfhcal_ecal1 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_proton/%1.1f/G4EICDetector_g4femc_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_proton_fecalfhcal_hcal2 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_proton/%1.1f/G4EICDetector_g4fhcal_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      TTree *const t_proton_fecalfhcal_ecal2 =  (TTree *) (new TFile(Form("/media/nschmidt/local/AnalysisSoftwareEIC/resolutionHCAL/EIC_simulation/FECAL_FHCAL_proton/%1.1f/G4EICDetector_g4femc_eval.root",partE[epart]), "READ"))->Get("ntp_gshower");
      if(t_proton_fecalfhcal_hcal1){
        activeEboth_proton++;
        h_proton_ecalfhcal_E[epart] 	= new TH1F(Form("h_proton_ecalfhcal_E_%d",epart), "", 200, 0.0, 2.0);

        fill_histogram(h_proton_ecalfhcal_E[epart], t_proton_fecalfhcal_hcal1, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_proton_fecalfhcal_ecal1){
        h_proton_ecal_E[epart] 	= new TH1F(Form("h_proton_ecal_E_%d",epart), "", 200, 0.0, 2.0);

        fill_histogram(h_proton_ecal_E[epart], t_proton_fecalfhcal_ecal1, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_proton_fecalfhcal_hcal2 && t_proton_fecalfhcal_ecal2){
        h2D_NTowers_HCAL_proton[epart] 	= new TH2F(Form("h2D_NTowers_HCAL_proton%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_NTowers_ECAL_proton[epart] 	= new TH2F(Form("h2D_NTowers_ECAL_proton%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_proton_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_proton_ecalfhcal_E_%d",epart), "", 100, 0,partE[epart]*1.15, 100, 0,partE[epart]*1.15);
        h_proton_ecalfhcalcomb_E[epart] 	= new TH1F(Form("h_proton_ecalfhcalcomb_E_%d",epart), "", 200, 0.0, 2.0);
        h_proton_ecalfhcalcomb_E_calibrated[epart] 	= new TH1F(Form("h_proton_ecalfhcalcomb_E_calibrated_%d",epart), "", 200, 0.0, 2.0);
        Float_t measured_energy1,measured_energy2;
        t_proton_fecalfhcal_hcal2->SetBranchAddress("e", &measured_energy1);
        t_proton_fecalfhcal_ecal2->SetBranchAddress("e", &measured_energy2);
        Float_t true_eta1,true_eta2;
        t_proton_fecalfhcal_hcal2->SetBranchAddress("geta", &true_eta1);
        t_proton_fecalfhcal_ecal2->SetBranchAddress("geta", &true_eta2);
        Float_t nTowers1,nTowers2;
        t_proton_fecalfhcal_hcal2->SetBranchAddress("ntowers", &nTowers1);
        t_proton_fecalfhcal_ecal2->SetBranchAddress("ntowers", &nTowers2);

        Int_t nentries = Int_t(t_proton_fecalfhcal_hcal2->GetEntries()) > Int_t(t_proton_fecalfhcal_ecal2->GetEntries()) ? Int_t(t_proton_fecalfhcal_ecal2->GetEntries()) : Int_t(t_proton_fecalfhcal_hcal2->GetEntries());
        Float_t lastEvt = 0;
        for (Int_t i = 0; i < nentries; ++i) {
          if (t_proton_fecalfhcal_hcal2->LoadTree(i) < 0 || t_proton_fecalfhcal_ecal2->LoadTree(i) < 0)
            break;

          t_proton_fecalfhcal_hcal2->GetEntry(i);
          t_proton_fecalfhcal_ecal2->GetEntry(i);
          if (((true_eta1 > 1.7 && true_eta1 < 3.0)
            || (true_eta2 > 1.7 && true_eta2 < 3.0))  //1.7 < 3
          && (measured_energy1 > 0.35) // above MIP
          && (measured_energy1 > 0.35) // above MIP
          && nTowers1>=2 && nTowers2>=2 // actual cluster with >1 tower
          ){
            h2D_NTowers_HCAL_proton[epart]->Fill(measured_energy1, nTowers1);
            h2D_NTowers_ECAL_proton[epart]->Fill(measured_energy2, nTowers2);
            h2D_proton_ecalfhcal_E[epart]->Fill(measured_energy2, measured_energy1);
            h_proton_ecalfhcalcomb_E[epart]->Fill((measured_energy2+measured_energy1)/partE[epart]);
            h_proton_ecalfhcalcomb_E_calibrated[epart]->Fill(((measured_energy2/calibconstECAL)+(measured_energy1/FunctionHCAL_calib->Eval(measured_energy1)))/partE[epart]);
          }
        }
        h_proton_ecalfhcalcomb_E[epart]->Scale(1 / h_proton_ecalfhcalcomb_E[epart]->GetEntries());
        h_proton_ecalfhcalcomb_E_calibrated[epart]->Scale(1 / h_proton_ecalfhcalcomb_E_calibrated[epart]->GetEntries());
      }
    }
  }

  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
	split_canvas(cPNG, "cPNG1", activeE);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 1.5,1000,0, 0.070);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_reso[nEne] = {NULL};
  TF1 *fit_reso_NCell[nEne] = {NULL};
  TF1 *fit_reso_eta1[nEne] = {NULL};
  TF1 *fit_reso_eta2[nEne] = {NULL};
  TF1 *fit_reso_eta3[nEne] = {NULL};
  TF1 *fit_reso_eta4[nEne] = {NULL};
  TF1 *fit_reso_gamma[nEne] = {NULL};

  TH1F*  h_pion_resolution_fhcal_E 	= new TH1F("h_pion_resolution_fhcal_E", "", 400, 0.0, 200.0);
  TH1F*  h_pion_resolution_fhcal_1oE 	= new TH1F("h_pion_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  Int_t padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummy->DrawCopy();
    DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
    h_pion_fhcal_E[epart]->Draw("same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    fit_reso[epart] = new TF1(Form("fit_reso_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);

    h_pion_fhcal_E[epart]->Fit(fit_reso[epart],"L0RMEQ");
    fit_reso[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
    fit_reso[epart]->Draw("same");

    drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    h_pion_resolution_fhcal_E->SetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(2)/fit_reso[epart]->GetParameter(2),2)+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2))); //+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)
    h_pion_resolution_fhcal_1oE->SetBinContent(h_pion_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_resolution.%s", outputDir.Data(), suffix.Data()));
  
  TH2F * histEPlotDummyGammaPion                           = new TH2F("histEPlotDummyGammaPion","histEPlotDummyGammaPion",1000,0, 1.5,1000,0, 0.13);
  SetStyleHistoTH2ForGraphs(histEPlotDummyGammaPion, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  histEPlotDummyGammaPion->GetXaxis()->SetNoExponent();
  histEPlotDummyGammaPion->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyGammaPion->GetXaxis()->SetMoreLogLabels(kTRUE);

  split_canvas(cPNG, "cPNG2", activeEgamma);

  TH1F*  h_gamma_resolution_fhcal_E 	= new TH1F("h_gamma_resolution_fhcal_E", "", 400, 0.0, 200.0);
  TH1F*  h_gamma_resolution_fhcal_1oE 	= new TH1F("h_gamma_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  TH1F* h_gamma_hcal_mean_E 	    = new TH1F("h_gamma_hcal_mean_E", "", 400, 0.0, 200.0);
  TH1F* h_gamma_ecal_mean_E 	    = new TH1F("h_gamma_ecal_mean_E", "", 400, 0.0, 200.0);
  TH1F* h_gamma_ecal_mean_E_calibrated 	    = new TH1F("h_gamma_ecal_mean_E_calibrated", "", 400, 0.0, 200.0);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart] && !useEgamma[epart]) continue;
    if(!useEgamma[epart]) continue;
    if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyGammaPion->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV}",partE[epart]),0.90,0.22,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("HCAL",0.90,0.16,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart]){
      DrawGammaSetMarker(h_pion_fhcal_E_allclus[epart], 24, 2, kGray+1, kGray+1);
      h_pion_fhcal_E_allclus[epart]->Draw("same");
      DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      h_pion_fhcal_E[epart]->Draw("same");
      DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      fit_reso[epart]->Draw("same");
      // fit_reso[epart]->SetRange(0.6,1.2);
      // DrawGammaSetMarkerTF1(fit_reso[epart], 2, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd("#color[602]{#pi^{-}}",0.25,0.60,2.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    if(useEgamma[epart]){
      DrawGammaSetMarker(h_gamma_fhcal_E[epart], 20, 2, kGreen+2, kGreen+2);
      h_gamma_fhcal_E[epart]->Draw("same");

      fit_reso_gamma[epart] = new TF1(Form("fit_reso_%d",epart), "gaus", 0.75,1.0);
      h_gamma_fhcal_E[epart]->Fit(fit_reso_gamma[epart],"L0RMEQ");
      fit_reso_gamma[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_gamma[epart], 1, 2, kBlack);
      fit_reso_gamma[epart]->Draw("same");

      drawLatexAdd(Form("mean: %1.2f ",fit_reso_gamma[epart]->GetParameter(1)),0.9,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gamma[epart]->GetParameter(2)),0.9,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gamma[epart]->GetParameter(2)/fit_reso_gamma[epart]->GetParameter(1)),0.9,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#color[418]{#gamma}",0.85,0.60,2.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      h_gamma_hcal_mean_E->SetBinContent(h_gamma_hcal_mean_E->FindBin(partE[epart]*fit_reso_gamma[epart]->GetParameter(1)),fit_reso_gamma[epart]->GetParameter(1));
      // h_gamma_hcal_mean_E->SetBinContent(h_gamma_hcal_mean_E->FindBin(partE[epart]),fit_reso_gamma[epart]->GetParameter(1));
      h_gamma_resolution_fhcal_E->SetBinContent(h_gamma_resolution_fhcal_E->FindBin(partE[epart]),100*fit_reso_gamma[epart]->GetParameter(2)/fit_reso_gamma[epart]->GetParameter(1));
      h_gamma_resolution_fhcal_1oE->SetBinContent(h_gamma_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso_gamma[epart]->GetParameter(2)/fit_reso_gamma[epart]->GetParameter(1));
    }
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_gamma_pion_resolution.%s", outputDir.Data(), suffix.Data()));
  
  split_canvas(cPNG, "cPNG3", activeE);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart] ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyGammaPion->GetYaxis()->SetRangeUser(0.,0.099);
    histEPlotDummyGammaPion->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV}",partE[epart]),0.90,0.22,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{-} on HCAL",0.90,0.16,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart]){
      DrawGammaSetMarker(h_pion_fhcal_E_allclus[epart], 24, 2, kGray+1, kGray+1);
      h_pion_fhcal_E_allclus[epart]->Draw("same");

      fit_reso_NCell[epart] = new TF1(Form("fit_reso_NCell_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);

      h_pion_fhcal_E_allclus[epart]->Fit(fit_reso_NCell[epart],"L0RMEQ");
      fit_reso_NCell[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_NCell[epart], 1, 2, kGray+3);
      fit_reso_NCell[epart]->Draw("same");

      drawLatexAdd(Form("mean: %1.2f ",fit_reso_NCell[epart]->GetParameter(1)),0.9,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso_NCell[epart]->GetParameter(2)),0.9,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_NCell[epart]->GetParameter(2)/fit_reso_NCell[epart]->GetParameter(1)),0.9,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

      DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      h_pion_fhcal_E[epart]->Draw("same");
      DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      fit_reso[epart]->Draw("same");
      // fit_reso[epart]->SetRange(0.6,1.2);
      // DrawGammaSetMarkerTF1(fit_reso[epart], 2, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd("#color[923]{all clusters}",0.9,0.60,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#color[602]{#it{N}_{cell} > 1}",0.20,0.60,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_pion_NCell_resolution.%s", outputDir.Data(), suffix.Data()));
  
  
  split_canvas(cPNG, "cPNG3", activeE);

  TH1F*  h_cpion_eta1_resolution_fhcal_1oE 	= new TH1F("h_cpion_eta1_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  TH1F*  h_cpion_eta2_resolution_fhcal_1oE 	= new TH1F("h_cpion_eta2_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  TH1F*  h_cpion_eta3_resolution_fhcal_1oE 	= new TH1F("h_cpion_eta3_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  TH1F*  h_cpion_eta4_resolution_fhcal_1oE 	= new TH1F("h_cpion_eta4_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart] ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyGammaPion->GetYaxis()->SetRangeUser(0.,0.099);
    histEPlotDummyGammaPion->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV}",partE[epart]),0.90,0.22,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{-} on HCAL",0.90,0.16,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart]){
      DrawGammaSetMarker(h_pion_fhcal_E_allclus[epart], 24, 2, kGray+1, kGray+1);
      h_pion_fhcal_E_allclus[epart]->Draw("same");
      fit_reso_NCell[epart] = new TF1(Form("fit_reso_NCell_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);
      h_pion_fhcal_E_allclus[epart]->Fit(fit_reso_NCell[epart],"L0RMEQ");
      fit_reso_NCell[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_NCell[epart], 1, 2, kGray+3);
      fit_reso_NCell[epart]->Draw("same");

      DrawGammaSetMarker(h_pion_fhcal_E_eta1[epart], 24, 2, kCyan+2, kCyan+2);
      h_pion_fhcal_E_eta1[epart]->Draw("same");
      fit_reso_eta1[epart] = new TF1(Form("fit_reso_eta1_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);
      h_pion_fhcal_E_eta1[epart]->Fit(fit_reso_eta1[epart],"L0RMEQ");
      fit_reso_eta1[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_eta1[epart], 1, 2, kCyan+4);
      fit_reso_eta1[epart]->Draw("same");

      DrawGammaSetMarker(h_pion_fhcal_E_eta2[epart], 24, 2, kBlue+2, kBlue+2);
      h_pion_fhcal_E_eta2[epart]->Draw("same");
      fit_reso_eta2[epart] = new TF1(Form("fit_reso_eta2_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);
      h_pion_fhcal_E_eta2[epart]->Fit(fit_reso_eta2[epart],"L0RMEQ");
      fit_reso_eta2[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_eta2[epart], 1, 2, kBlue+4);
      fit_reso_eta2[epart]->Draw("same");

      DrawGammaSetMarker(h_pion_fhcal_E_eta3[epart], 24, 2, kMagenta+2, kMagenta+2);
      h_pion_fhcal_E_eta3[epart]->Draw("same");
      fit_reso_eta3[epart] = new TF1(Form("fit_reso_eta3_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);
      h_pion_fhcal_E_eta3[epart]->Fit(fit_reso_eta3[epart],"L0RMEQ");
      fit_reso_eta3[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_eta3[epart], 1, 2, kMagenta+4);
      fit_reso_eta3[epart]->Draw("same");

      DrawGammaSetMarker(h_pion_fhcal_E_eta4[epart], 24, 2, kRed+2, kRed+2);
      h_pion_fhcal_E_eta4[epart]->Draw("same");
      fit_reso_eta4[epart] = new TF1(Form("fit_reso_eta4_%d",epart), "gaus", fitrange_low[epart]*0.8,fitrange_hig[epart]);
      h_pion_fhcal_E_eta4[epart]->Fit(fit_reso_eta4[epart],"L0RMEQ");
      fit_reso_eta4[epart]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_eta4[epart], 1, 2, kRed+4);
      fit_reso_eta4[epart]->Draw("same");

      h_cpion_eta1_resolution_fhcal_1oE->SetBinContent(h_cpion_eta1_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso_eta1[epart]->GetParameter(2)/fit_reso_eta1[epart]->GetParameter(1));
      h_cpion_eta2_resolution_fhcal_1oE->SetBinContent(h_cpion_eta2_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso_eta2[epart]->GetParameter(2)/fit_reso_eta2[epart]->GetParameter(1));
      h_cpion_eta3_resolution_fhcal_1oE->SetBinContent(h_cpion_eta3_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso_eta3[epart]->GetParameter(2)/fit_reso_eta3[epart]->GetParameter(1));
      h_cpion_eta4_resolution_fhcal_1oE->SetBinContent(h_cpion_eta4_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso_eta4[epart]->GetParameter(2)/fit_reso_eta4[epart]->GetParameter(1));


  TLegend* legendEtaCpion           = GetAndSetLegend2(0.45, 0.93-(6*textSizeLabelsRel), 0.9, 0.93,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendEtaCpion->AddEntry(h_pion_fhcal_E_allclus[epart],Form("1.4<#eta<3.5: #sigma/#it{E}: %1.2f %%",100*fit_reso_NCell[epart]->GetParameter(2)/fit_reso_NCell[epart]->GetParameter(1)),"pl");
  legendEtaCpion->AddEntry(h_pion_fhcal_E_eta1[epart],Form("1.2<#eta<1.7: #sigma/#it{E}: %1.2f %%",100*fit_reso_eta1[epart]->GetParameter(2)/fit_reso_eta1[epart]->GetParameter(1)),"pl");
  legendEtaCpion->AddEntry(h_pion_fhcal_E_eta2[epart],Form("1.7<#eta<2.5: #sigma/#it{E}: %1.2f %%",100*fit_reso_eta2[epart]->GetParameter(2)/fit_reso_eta2[epart]->GetParameter(1)),"pl");
  legendEtaCpion->AddEntry(h_pion_fhcal_E_eta3[epart],Form("2.5<#eta<3.0: #sigma/#it{E}: %1.2f %%",100*fit_reso_eta3[epart]->GetParameter(2)/fit_reso_eta3[epart]->GetParameter(1)),"pl");
  legendEtaCpion->AddEntry(h_pion_fhcal_E_eta4[epart],Form("3.0<#eta<4.0: #sigma/#it{E}: %1.2f %%",100*fit_reso_eta4[epart]->GetParameter(2)/fit_reso_eta4[epart]->GetParameter(1)),"pl");
  legendEtaCpion->Draw();

      // DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      // h_pion_fhcal_E[epart]->Draw("same");
      // DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      // // fit_reso[epart]->SetRange(0.6,1.2);
      // // DrawGammaSetMarkerTF1(fit_reso[epart], 2, 2, kBlack);
      // // fit_reso[epart]->Draw("same");
      // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd("#color[923]{all clusters}",0.9,0.60,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      // drawLatexAdd("#color[602]{#it{N}_{cell} > 1}",0.20,0.60,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_pion_eta_resolution.%s", outputDir.Data(), suffix.Data()));
  


  TCanvas* cPNGBoth;
  split_canvas(cPNGBoth, "cPNGBoth", activeEboth);

  TH2F * histEPlotDummyBoth                           = new TH2F("histEPlotDummyBoth","histEPlotDummyBoth",1000,0, 1.5,1000,0, 0.085);
  SetStyleHistoTH2ForGraphs(histEPlotDummyBoth, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummyBoth->GetXaxis()->SetNoExponent();
  histEPlotDummyBoth->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyBoth->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_resoBothComb[nEne] = {NULL};
  TF1 *fit_resoBothComb_calib[nEne] = {NULL};
  TF1 *fit_resoBothComb_calibrated[nEne] = {NULL};
  TF1 *fit_resoBoth[nEne] = {NULL};
  TF1 *fit_resoBothEMC[nEne] = {NULL};
  TLegend* legendHE_pion = NULL;
  TH1F*  h_pion_resolution_fecalfhcal_E 	= new TH1F("h_pion_resolution_fecalfhcal_E", "", 200, 0.0, 100.0);
  TH1F*  h_pion_resolution_fecalfhcal_1oE 	= new TH1F("h_pion_resolution_fecalfhcal_1oE", "", 500, 0.0, 1.0);
  TH1F*  h_pion_resolution_comb_E 	= new TH1F("h_pion_resolution_comb_E", "", 200, 0.0, 100.0);
  TH1F*  h_pion_resolution_comb_1oE 	= new TH1F("h_pion_resolution_comb_1oE", "", 500, 0.0, 1.0);

  TH1F* h_cpion_comb_mean_E 	    = new TH1F("h_cpion_comb_mean_E", "", 400, 0.0, 200.0);
  TH1F* h_cpion_comb_mean_E_calib 	    = new TH1F("h_cpion_comb_mean_E_calib", "", 400, 0.0, 200.0);
  TH1F* h_cpion_comb_mean_E_calibrated 	    = new TH1F("h_cpion_comb_mean_E_calibrated", "", 400, 0.0, 200.0);
  TH1F* h_proton_comb_mean_E 	    = new TH1F("h_proton_comb_mean_E", "", 400, 0.0, 200.0);
  TH1F* h_proton_comb_mean_E_calibrated 	    = new TH1F("h_proton_comb_mean_E_calibrated", "", 400, 0.0, 200.0);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    if(!h_pion_ecal_E[epart]) continue;
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    DrawGammaSetMarker(h_pion_ecal_E[epart], 28, 2, kGreen+2, kGreen+2);
    h_pion_ecal_E[epart]->Draw("same");
    fit_resoBothEMC[epart] = new TF1(Form("fit_resoBothEMC_%d",epart), "gaus", 0.2,0.6);
    h_pion_ecal_E[epart]->Fit(fit_resoBothEMC[epart],"L0RMEQ");
    fit_resoBothEMC[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_resoBothEMC[epart], 1, 2, kGreen+4);
    fit_resoBothEMC[epart]->Draw("same");

    DrawGammaSetMarker(h_pion_ecalfhcal_E[epart], 27, 2, kBlue+2, kBlue+2);
    h_pion_ecalfhcal_E[epart]->Draw("same");
    fit_resoBoth[epart] = new TF1(Form("fit_resoBoth_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);
    h_pion_ecalfhcal_E[epart]->Fit(fit_resoBoth[epart],"L0RMEQ");
    fit_resoBoth[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_resoBoth[epart], 1, 2, kBlack);
    fit_resoBoth[epart]->Draw("same");

    DrawGammaSetMarker(h_pion_ecalfhcalcomb_E[epart], 20, 2, kRed+2, kRed+2);
    h_pion_ecalfhcalcomb_E[epart]->Draw("same");
    fit_resoBothComb[epart] = new TF1(Form("fit_resoBothComb_%d",epart), "gaus", partE[epart]<50 ? 0.35+0.008*partE[epart] : 0.65,0.9);
    h_pion_ecalfhcalcomb_E[epart]->Fit(fit_resoBothComb[epart],"L0RMEQ");
    fit_resoBothComb[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_resoBothComb[epart], 1, 2,  kRed+4);
    fit_resoBothComb[epart]->Draw("same");

    DrawGammaSetMarker(h_pion_ecalfhcalcomb_E_calib[epart], 20, 2, kMagenta+2, kMagenta+2);
    // h_pion_ecalfhcalcomb_E_calib[epart]->Draw("same");
    fit_resoBothComb_calib[epart] = new TF1(Form("fit_resoBothComb_calib%d",epart), "gaus", partE[epart]<50 ? 0.3+0.008*partE[epart] : 0.65,0.9);
    h_pion_ecalfhcalcomb_E_calib[epart]->Fit(fit_resoBothComb_calib[epart],"L0RMEQ");
    fit_resoBothComb_calib[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_resoBothComb_calib[epart], 1, 2,  kMagenta+4);
    // fit_resoBothComb_calib[epart]->Draw("same");

    DrawGammaSetMarker(h_pion_ecalfhcalcomb_E_calibrated[epart], 20, 2, kOrange+2, kOrange+2);
    h_pion_ecalfhcalcomb_E_calibrated[epart]->Draw("same");
    fit_resoBothComb_calibrated[epart] = new TF1(Form("fit_resoBothComb_calibrated%d",epart), "gaus", partE[epart]<50 ? 0.4+0.008*partE[epart] : 0.85,partE[epart]<30 ? 1.3 : 1.15);
    h_pion_ecalfhcalcomb_E_calibrated[epart]->Fit(fit_resoBothComb_calibrated[epart],"L0RMEQ");
    fit_resoBothComb_calibrated[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_resoBothComb_calibrated[epart], 1, 2,  kOrange+4);
    fit_resoBothComb_calibrated[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV #pi^{-}",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    legendHE_pion           = GetAndSetLegend2(0.20, 0.92-(5*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendHE_pion->AddEntry(h_pion_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBothEMC[epart]->GetParameter(1),fit_resoBothEMC[epart]->GetParameter(2),100*fit_resoBothEMC[epart]->GetParameter(2)/fit_resoBothEMC[epart]->GetParameter(1)),"pl");
    legendHE_pion->AddEntry(h_pion_ecalfhcal_E[epart],Form("HCAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBoth[epart]->GetParameter(1),fit_resoBoth[epart]->GetParameter(2),100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1)),"pl");
    legendHE_pion->AddEntry(h_pion_ecalfhcalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBothComb[epart]->GetParameter(1),fit_resoBothComb[epart]->GetParameter(2),100*fit_resoBothComb[epart]->GetParameter(2)/fit_resoBothComb[epart]->GetParameter(1)),"pl");
    legendHE_pion->AddEntry(h_pion_ecalfhcalcomb_E_calibrated[epart],Form("Calib: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBothComb_calibrated[epart]->GetParameter(1),fit_resoBothComb_calibrated[epart]->GetParameter(2),100*fit_resoBothComb_calibrated[epart]->GetParameter(2)/fit_resoBothComb_calibrated[epart]->GetParameter(1)),"pl");
    // legendNonLin->AddEntry((TObject*)0,"","");
    legendHE_pion->Draw();

    // drawLatexAdd(Form("mean: %1.2f ",fit_resoBoth[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("sigma: %1.2f ",fit_resoBoth[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // h_cpion_comb_mean_E->SetBinContent(h_cpion_comb_mean_E->FindBin(partE[epart]),fit_resoBoth[epart]->GetParameter(1));
    // h_cpion_comb_mean_E_calib->SetBinContent(h_cpion_comb_mean_E_calib->FindBin(partE[epart]),fit_resoBothComb_calib[epart]->GetParameter(1));
    // h_cpion_comb_mean_E_calibrated->SetBinContent(h_cpion_comb_mean_E_calibrated->FindBin(partE[epart]),fit_resoBothComb_calibrated[epart]->GetParameter(1));
    h_cpion_comb_mean_E->SetBinContent(h_cpion_comb_mean_E->FindBin(partE[epart]*fit_resoBoth[epart]->GetParameter(1)),fit_resoBoth[epart]->GetParameter(1));
    h_cpion_comb_mean_E_calib->SetBinContent(h_cpion_comb_mean_E_calib->FindBin(partE[epart]*fit_resoBothComb_calib[epart]->GetParameter(1)),fit_resoBothComb_calib[epart]->GetParameter(1));
    h_cpion_comb_mean_E_calibrated->SetBinContent(h_cpion_comb_mean_E_calibrated->FindBin(partE[epart]*fit_resoBothComb_calibrated[epart]->GetParameter(1)),fit_resoBothComb_calibrated[epart]->GetParameter(1));
    h_pion_resolution_fecalfhcal_E->SetBinContent(h_pion_resolution_fecalfhcal_E->FindBin(partE[epart]),100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1));
    // h_pion_resolution_fecalfhcal_E->SetBinError(h_pion_resolution_fecalfhcal_E->FindBin(partE[epart]),100*TMath::Sqrt(TMath::Power(fit_resoBoth[epart]->GetParError(2)/fit_resoBoth[epart]->GetParameter(2),2)+TMath::Power(fit_resoBoth[epart]->GetParError(1)/fit_resoBoth[epart]->GetParameter(1),2)));
    h_pion_resolution_fecalfhcal_1oE->SetBinContent(h_pion_resolution_fecalfhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1));
    // h_pion_resolution_fecalfhcal_1oE->SetBinError(h_pion_resolution_fecalfhcal_1oE->FindBin(partE[epart]),100*TMath::Sqrt(TMath::Power(fit_resoBoth[epart]->GetParError(2)/fit_resoBoth[epart]->GetParameter(2),2)+TMath::Power(fit_resoBoth[epart]->GetParError(1)/fit_resoBoth[epart]->GetParameter(1),2)));
    h_pion_resolution_comb_E->SetBinContent(h_pion_resolution_comb_E->FindBin(partE[epart]),100*fit_resoBothComb[epart]->GetParameter(2)/fit_resoBothComb[epart]->GetParameter(1));
    h_pion_resolution_comb_1oE->SetBinContent(h_pion_resolution_comb_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_resoBothComb[epart]->GetParameter(2)/fit_resoBothComb[epart]->GetParameter(1));
    padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_HCAL_resolution_cpion.%s", outputDir.Data(), suffix.Data()));

  if(energySingleCalibPlot!=-1){
   split_canvas(cPNGBoth, "cPNGBoth", 1);
   padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    if(!h_pion_ecal_E[epart]) continue;
    if( energySingleCalibPlot!=partE[epart]) continue;
  	cPNGBoth->SetLeftMargin(0.1);
  	cPNGBoth->SetBottomMargin(0.1);
  	cPNGBoth->SetTopMargin(0.01);
  	cPNGBoth->SetRightMargin(0.01);
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    h_pion_ecal_E[epart]->Draw("same");
    fit_resoBothEMC[epart]->Draw("same");

    h_pion_ecalfhcal_E[epart]->Draw("same");
    fit_resoBoth[epart]->Draw("same");

    h_pion_ecalfhcalcomb_E[epart]->Draw("same");
    fit_resoBothComb[epart]->Draw("same");

    h_pion_ecalfhcalcomb_E_calib[epart]->Draw("same");
    fit_resoBothComb_calib[epart]->Draw("same");

    h_pion_ecalfhcalcomb_E_calibrated[epart]->Draw("same");
    fit_resoBothComb_calibrated[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV #pi^{-}",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    legendHE_pion           = GetAndSetLegend2(0.20, 0.92-(5*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendHE_pion->AddEntry(h_pion_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBothEMC[epart]->GetParameter(1),fit_resoBothEMC[epart]->GetParameter(2),100*fit_resoBothEMC[epart]->GetParameter(2)/fit_resoBothEMC[epart]->GetParameter(1)),"pl");
    legendHE_pion->AddEntry(h_pion_ecalfhcal_E[epart],Form("HCAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBoth[epart]->GetParameter(1),fit_resoBoth[epart]->GetParameter(2),100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1)),"pl");
    legendHE_pion->AddEntry(h_pion_ecalfhcalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBothComb[epart]->GetParameter(1),fit_resoBothComb[epart]->GetParameter(2),100*fit_resoBothComb[epart]->GetParameter(2)/fit_resoBothComb[epart]->GetParameter(1)),"pl");
    legendHE_pion->AddEntry(h_pion_ecalfhcalcomb_E_calibrated[epart],Form("Calib: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBothComb_calibrated[epart]->GetParameter(1),fit_resoBothComb_calibrated[epart]->GetParameter(2),100*fit_resoBothComb_calibrated[epart]->GetParameter(2)/fit_resoBothComb_calibrated[epart]->GetParameter(1)),"pl");
    legendHE_pion->Draw();
    padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_HCAL_resolution_cpion_single.%s", outputDir.Data(), suffix.Data()));
  }
  split_canvas(cPNGBoth, "cPNGBoth", activeEboth_gamma);
  histEPlotDummyBoth                           = new TH2F("histEPlotDummyBoth","histEPlotDummyBoth",1000,0, 1.5,1000,0, 0.265);
  SetStyleHistoTH2ForGraphs(histEPlotDummyBoth, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummyBoth->GetXaxis()->SetNoExponent();
  histEPlotDummyBoth->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyBoth->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_gamma_resoBothComb[nEne] = {NULL};
  TF1 *fit_gamma_resoBoth[nEne] = {NULL};
  TF1 *fit_gamma_resoBothEMC[nEne] = {NULL};
  TF1 *fit_gamma_resoBothEMC_calibrated[nEne] = {NULL};
  TLegend* legendHE_gamma = NULL;
  TH1F*  h_gamma_resolution_fecalfhcal_E 	= new TH1F("h_gamma_resolution_fecalfhcal_E", "", 200, 0.0, 100.0);
  TH1F*  h_gamma_resolution_fecalfhcal_1oE 	= new TH1F("h_gamma_resolution_fecalfhcal_1oE", "", 500, 0.0, 1.0);
  TH1F*  h_gamma_resolution_fecal_E 	= new TH1F("h_gamma_resolution_fecal_E", "", 200, 0.0, 100.0);
  TH1F*  h_gamma_resolution_fecal_1oE 	= new TH1F("h_gamma_resolution_fecal_1oE", "", 500, 0.0, 1.0);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothgamma[epart]) continue;
    if(!h_gamma_ecal_E[epart]) continue;
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    DrawGammaSetMarker(h_gamma_ecal_E[epart], 28, 2, kGreen+2, kGreen+2);
    h_gamma_ecal_E[epart]->Draw("same");
    fit_gamma_resoBothEMC[epart] = new TF1(Form("fit_gamma_resoBothEMC_%d",epart), "gaus", 0.8,0.9);
    h_gamma_ecal_E[epart]->Fit(fit_gamma_resoBothEMC[epart],"L0RMEQ");
    fit_gamma_resoBothEMC[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_gamma_resoBothEMC[epart], 1, 2, kGreen+4);
    fit_gamma_resoBothEMC[epart]->Draw("same");
    h_gamma_resolution_fecal_E->SetBinContent(h_gamma_resolution_fecal_E->FindBin(partE[epart]),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1));
    h_gamma_resolution_fecal_1oE->SetBinContent(h_gamma_resolution_fecal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1));
    DrawGammaSetMarker(h_gamma_ecal_E_calibrated[epart], 28, 2, kOrange+2, kOrange+2);
    h_gamma_ecal_E_calibrated[epart]->Draw("same");
    fit_gamma_resoBothEMC_calibrated[epart] = new TF1(Form("fit_gamma_resoBothEMC_calibrated%d",epart), "gaus", 0.95,1.05);
    h_gamma_ecal_E_calibrated[epart]->Fit(fit_gamma_resoBothEMC_calibrated[epart],"L0RMEQ");
    fit_gamma_resoBothEMC_calibrated[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_gamma_resoBothEMC_calibrated[epart], 1, 2, kOrange+4);
    fit_gamma_resoBothEMC_calibrated[epart]->Draw("same");
    h_gamma_ecal_mean_E->SetBinContent(h_gamma_ecal_mean_E->FindBin(partE[epart]*fit_gamma_resoBothEMC[epart]->GetParameter(1)),fit_gamma_resoBothEMC[epart]->GetParameter(1));
    // h_gamma_ecal_mean_E->SetBinContent(h_gamma_ecal_mean_E->FindBin(partE[epart]),fit_gamma_resoBothEMC[epart]->GetParameter(1));
    h_gamma_ecal_mean_E_calibrated->SetBinContent(h_gamma_ecal_mean_E_calibrated->FindBin(partE[epart]*fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(1)),fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(1));
    // h_gamma_ecal_mean_E_calibrated->SetBinContent(h_gamma_ecal_mean_E_calibrated->FindBin(partE[epart]),fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(1));

    DrawGammaSetMarker(h_gamma_ecalfhcal_E[epart], 27, 2, kBlue+2, kBlue+2);
    h_gamma_ecalfhcal_E[epart]->Draw("same");
    fit_gamma_resoBoth[epart] = new TF1(Form("fit_gamma_resoBoth_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);
    h_gamma_ecalfhcal_E[epart]->Fit(fit_gamma_resoBoth[epart],"L0RMEQ");
    fit_gamma_resoBoth[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_gamma_resoBoth[epart], 1, 2, kBlack);
    fit_gamma_resoBoth[epart]->Draw("same");

    DrawGammaSetMarker(h_gamma_ecalfhcalcomb_E[epart], 20, 2, kRed+2, kRed+2);
    h_gamma_ecalfhcalcomb_E[epart]->Draw("same");
    fit_gamma_resoBothComb[epart] = new TF1(Form("fit_gamma_resoBothComb_%d",epart), "gaus", partE[epart]==60 ? 0.83 : 0.8,partE[epart]==60 ? 0.87 : 0.9); //6+0.005*partE[epart]
    h_gamma_ecalfhcalcomb_E[epart]->Fit(fit_gamma_resoBothComb[epart],"L0RMEQ");
    fit_gamma_resoBothComb[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_gamma_resoBothComb[epart], 1, 2,  kRed+4);
    fit_gamma_resoBothComb[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV #gamma",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    legendHE_gamma           = GetAndSetLegend2(0.20, 0.92-(5*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendHE_gamma->AddEntry(h_gamma_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothEMC[epart]->GetParameter(1),fit_gamma_resoBothEMC[epart]->GetParameter(2),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1)),"pl");
    legendHE_gamma->AddEntry(h_gamma_ecalfhcal_E[epart],Form("HCAL: no fit"),"pl");
    // legendHE_gamma->AddEntry(h_gamma_ecalfhcal_E[epart],Form("HCAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBoth[epart]->GetParameter(1),fit_gamma_resoBoth[epart]->GetParameter(0),100*fit_gamma_resoBoth[epart]->GetParameter(0)/fit_gamma_resoBoth[epart]->GetParameter(1)),"pl");
    legendHE_gamma->AddEntry(h_gamma_ecalfhcalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothComb[epart]->GetParameter(1),fit_gamma_resoBothComb[epart]->GetParameter(2),100*fit_gamma_resoBothComb[epart]->GetParameter(2)/fit_gamma_resoBothComb[epart]->GetParameter(1)),"pl");
    legendHE_gamma->AddEntry(h_gamma_ecal_E_calibrated[epart],Form("Calib: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(1),fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(2),100*fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(2)/fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(1)),"pl");
    // legendNonLin->AddEntry((TObject*)0,"","");
    legendHE_gamma->Draw();

    // drawLatexAdd(Form("mean: %1.2f ",fit_gamma_resoBoth[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("sigma: %1.2f ",fit_gamma_resoBoth[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_gamma_resoBoth[epart]->GetParameter(2)/fit_gamma_resoBoth[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // h_pion_resolution_fecalfhcal_E->SetBinContent(h_pion_resolution_fecalfhcal_E->FindBin(partE[epart]),100*fit_gamma_resoBoth[epart]->GetParameter(2)/fit_gamma_resoBoth[epart]->GetParameter(1));
    // h_pion_resolution_fecalfhcal_1oE->SetBinContent(h_pion_resolution_fecalfhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_gamma_resoBoth[epart]->GetParameter(2)/fit_gamma_resoBoth[epart]->GetParameter(1));
    padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_HCAL_resolution_gamma.%s", outputDir.Data(), suffix.Data()));

  if(energySingleCalibPlot!=-1){
   split_canvas(cPNGBoth, "cPNGBothSinglegamm", 1);
   padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothgamma[epart]) continue;
    if(!h_gamma_ecal_E[epart]) continue;
    if(energySingleCalibPlot!=partE[epart]) continue;
  	cPNGBoth->SetLeftMargin(0.1);
  	cPNGBoth->SetBottomMargin(0.1);
  	cPNGBoth->SetTopMargin(0.01);
  	cPNGBoth->SetRightMargin(0.01);
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    h_gamma_ecal_E[epart]->Draw("same");
    fit_gamma_resoBothEMC[epart]->Draw("same");
    h_gamma_ecal_E_calibrated[epart]->Draw("same");
    fit_gamma_resoBothEMC_calibrated[epart]->Draw("same");
    h_gamma_ecalfhcal_E[epart]->Draw("same");
    fit_gamma_resoBoth[epart]->Draw("same");
    h_gamma_ecalfhcalcomb_E[epart]->Draw("same");
    fit_gamma_resoBothComb[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV #gamma",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    legendHE_gamma           = GetAndSetLegend2(0.20, 0.92-(5*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendHE_gamma->AddEntry(h_gamma_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothEMC[epart]->GetParameter(1),fit_gamma_resoBothEMC[epart]->GetParameter(2),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1)),"pl");
    legendHE_gamma->AddEntry(h_gamma_ecalfhcal_E[epart],Form("HCAL: no fit"),"pl");
    // legendHE_gamma->AddEntry(h_gamma_ecalfhcal_E[epart],Form("HCAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBoth[epart]->GetParameter(1),fit_gamma_resoBoth[epart]->GetParameter(0),100*fit_gamma_resoBoth[epart]->GetParameter(0)/fit_gamma_resoBoth[epart]->GetParameter(1)),"pl");
    legendHE_gamma->AddEntry(h_gamma_ecalfhcalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothComb[epart]->GetParameter(1),fit_gamma_resoBothComb[epart]->GetParameter(2),100*fit_gamma_resoBothComb[epart]->GetParameter(2)/fit_gamma_resoBothComb[epart]->GetParameter(1)),"pl");
    legendHE_gamma->AddEntry(h_gamma_ecal_E_calibrated[epart],Form("Calib: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(1),fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(2),100*fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(2)/fit_gamma_resoBothEMC_calibrated[epart]->GetParameter(1)),"pl");
    // legendNonLin->AddEntry((TObject*)0,"","");
    legendHE_gamma->Draw();
    padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_HCAL_resolution_gamma_single.%s", outputDir.Data(), suffix.Data()));
  }
  split_canvas(cPNGBoth, "cPNGBothProton", activeEboth_gamma);
  histEPlotDummyBoth                           = new TH2F("histEPlotDummyBoth","histEPlotDummyBoth",1000,0, 1.5,1000,0, 0.075);
  SetStyleHistoTH2ForGraphs(histEPlotDummyBoth, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummyBoth->GetXaxis()->SetNoExponent();
  histEPlotDummyBoth->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyBoth->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_proton_resoBothComb[nEne] = {NULL};
  TF1 *fit_proton_resoBothComb_calibrated[nEne] = {NULL};
  TF1 *fit_proton_resoBoth[nEne] = {NULL};
  TF1 *fit_proton_resoBothEMC[nEne] = {NULL};
  TLegend* legendHE_proton = NULL;
  TH1F*  h_proton_resolution_fecalfhcal_E 	= new TH1F("h_proton_resolution_fecalfhcal_E", "", 200, 0.0, 100.0);
  TH1F*  h_proton_resolution_fecalfhcal_1oE 	= new TH1F("h_proton_resolution_fecalfhcal_1oE", "", 500, 0.0, 1.0);
  TH1F*  h_proton_resolution_fecal_E 	= new TH1F("h_proton_resolution_fecal_E", "", 200, 0.0, 100.0);
  TH1F*  h_proton_resolution_fecal_1oE 	= new TH1F("h_proton_resolution_fecal_1oE", "", 500, 0.0, 1.0);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothproton[epart]) continue;
    if(!h_proton_ecal_E[epart]) continue;
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    DrawGammaSetMarker(h_proton_ecal_E[epart], 28, 2, kGreen+2, kGreen+2);
    h_proton_ecal_E[epart]->Draw("same");
    fit_proton_resoBothEMC[epart] = new TF1(Form("fit_proton_resoBothEMC_%d",epart), "gaus", 0.2,0.5);
    h_proton_ecal_E[epart]->Fit(fit_proton_resoBothEMC[epart],"L0RMEQ");
    fit_proton_resoBothEMC[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_proton_resoBothEMC[epart], 1, 2, kGreen+4);
    fit_proton_resoBothEMC[epart]->Draw("same");
    h_proton_resolution_fecal_E->SetBinContent(h_proton_resolution_fecal_E->FindBin(partE[epart]),100*fit_proton_resoBothEMC[epart]->GetParameter(2)/fit_proton_resoBothEMC[epart]->GetParameter(1));
    h_proton_resolution_fecal_1oE->SetBinContent(h_proton_resolution_fecal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_proton_resoBothEMC[epart]->GetParameter(2)/fit_proton_resoBothEMC[epart]->GetParameter(1));


    DrawGammaSetMarker(h_proton_ecalfhcal_E[epart], 27, 2, kBlue+2, kBlue+2);
    h_proton_ecalfhcal_E[epart]->Draw("same");
    fit_proton_resoBoth[epart] = new TF1(Form("fit_proton_resoBoth_%d",epart), "gaus", 0.3,0.8);
    h_proton_ecalfhcal_E[epart]->Fit(fit_proton_resoBoth[epart],"L0RMEQ");
    fit_proton_resoBoth[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_proton_resoBoth[epart], 1, 2, kBlack);
    fit_proton_resoBoth[epart]->Draw("same");

    DrawGammaSetMarker(h_proton_ecalfhcalcomb_E[epart], 20, 2, kRed+2, kRed+2);
    h_proton_ecalfhcalcomb_E[epart]->Draw("same");
    fit_proton_resoBothComb[epart] = new TF1(Form("fit_proton_resoBothComb_%d",epart), "gaus", partE[epart]<40 ? 0.2+0.005*partE[epart] : 0.6,partE[epart]<40 ?0.55+0.01*partE[epart] : 0.8); //6+0.005*partE[epart]
    h_proton_ecalfhcalcomb_E[epart]->Fit(fit_proton_resoBothComb[epart],"L0RMEQ");
    fit_proton_resoBothComb[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_proton_resoBothComb[epart], 1, 2,  kRed+4);
    fit_proton_resoBothComb[epart]->Draw("same");

    DrawGammaSetMarker(h_proton_ecalfhcalcomb_E_calibrated[epart], 20, 2, kOrange+2, kOrange+2);
    h_proton_ecalfhcalcomb_E_calibrated[epart]->Draw("same");
    fit_proton_resoBothComb_calibrated[epart] = new TF1(Form("fit_proton_resoBothComb_calibrated_%d",epart), "gaus", 0.4+0.005*partE[epart],1.05+0.01*partE[epart]); //6+0.005*partE[epart]
    h_proton_ecalfhcalcomb_E_calibrated[epart]->Fit(fit_proton_resoBothComb_calibrated[epart],"L0RMEQ");
    fit_proton_resoBothComb_calibrated[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_proton_resoBothComb_calibrated[epart], 1, 2,  kOrange+4);
    fit_proton_resoBothComb_calibrated[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV proton",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    legendHE_proton           = GetAndSetLegend2(0.20, 0.92-(5*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    legendHE_proton->AddEntry(h_proton_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_proton_resoBothEMC[epart]->GetParameter(1),fit_proton_resoBothEMC[epart]->GetParameter(2),100*fit_proton_resoBothEMC[epart]->GetParameter(2)/fit_proton_resoBothEMC[epart]->GetParameter(1)),"pl");
    // legendHE_proton->AddEntry(h_proton_ecalfhcal_E[epart],Form("HCAL: no fit"),"pl");
    legendHE_proton->AddEntry(h_proton_ecalfhcal_E[epart],Form("HCAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_proton_resoBoth[epart]->GetParameter(1),fit_proton_resoBoth[epart]->GetParameter(0),100*fit_proton_resoBoth[epart]->GetParameter(0)/fit_proton_resoBoth[epart]->GetParameter(1)),"pl");
    legendHE_proton->AddEntry(h_proton_ecalfhcalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_proton_resoBothComb[epart]->GetParameter(1),fit_proton_resoBothComb[epart]->GetParameter(2),100*fit_proton_resoBothComb[epart]->GetParameter(2)/fit_proton_resoBothComb[epart]->GetParameter(1)),"pl");
    legendHE_proton->AddEntry(h_proton_ecalfhcalcomb_E_calibrated[epart],Form("Calib: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_proton_resoBothComb_calibrated[epart]->GetParameter(1),fit_proton_resoBothComb_calibrated[epart]->GetParameter(2),100*fit_proton_resoBothComb_calibrated[epart]->GetParameter(2)/fit_proton_resoBothComb_calibrated[epart]->GetParameter(1)),"pl");
    // legendNonLin->AddEntry((TObject*)0,"","");
    legendHE_proton->Draw();

    // drawLatexAdd(Form("mean: %1.2f ",fit_proton_resoBoth[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("sigma: %1.2f ",fit_proton_resoBoth[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_proton_resoBoth[epart]->GetParameter(2)/fit_proton_resoBoth[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // h_pion_resolution_fecalfhcal_E->SetBinContent(h_pion_resolution_fecalfhcal_E->FindBin(partE[epart]),100*fit_proton_resoBoth[epart]->GetParameter(2)/fit_proton_resoBoth[epart]->GetParameter(1));
    h_proton_resolution_fecalfhcal_1oE->SetBinContent(h_proton_resolution_fecalfhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_proton_resoBothComb[epart]->GetParameter(2)/fit_proton_resoBothComb[epart]->GetParameter(1));
    // h_proton_comb_mean_E->SetBinContent(h_proton_comb_mean_E->FindBin(partE[epart]),fit_proton_resoBothComb[epart]->GetParameter(1));
    // h_proton_comb_mean_E_calibrated->SetBinContent(h_proton_comb_mean_E_calibrated->FindBin(partE[epart]),fit_proton_resoBothComb_calibrated[epart]->GetParameter(1));
    h_proton_comb_mean_E->SetBinContent(h_proton_comb_mean_E->FindBin(partE[epart]*fit_proton_resoBothComb[epart]->GetParameter(1)),fit_proton_resoBothComb[epart]->GetParameter(1));
    h_proton_comb_mean_E_calibrated->SetBinContent(h_proton_comb_mean_E_calibrated->FindBin(partE[epart]*fit_proton_resoBothComb_calibrated[epart]->GetParameter(1)),fit_proton_resoBothComb_calibrated[epart]->GetParameter(1));
   padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_HCAL_resolution_proton.%s", outputDir.Data(), suffix.Data()));

  split_canvas(cPNGBoth, "cPNGBothProton", 1);
  histEPlotDummyBoth                           = new TH2F("histEPlotDummyBoth","histEPlotDummyBoth",1000,0.3, 1.05,1000,0, 0.145);
  SetStyleHistoTH2ForGraphs(histEPlotDummyBoth, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.1);
  histEPlotDummyBoth->GetXaxis()->SetNoExponent();
  histEPlotDummyBoth->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyBoth->GetXaxis()->SetMoreLogLabels(kTRUE);

  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothproton[epart]) continue;
    if(!h_proton_ecal_E[epart]) continue;
    if(partE[epart]!=60) continue;

    TCanvas* canvasECCE = new TCanvas("canvasECCE","",0,0,1100,1000);
    DrawGammaCanvasSettings( canvasECCE, 0.11, 0.01, 0.01, 0.105);
    histEPlotDummyBoth->DrawCopy();

    DrawGammaSetMarker(h_proton_ecalfhcalcomb_E[epart], 43, 4, kGray+2, kGray+2);
    h_proton_ecalfhcalcomb_E[epart]->Draw("same");
    DrawGammaSetMarkerTF1(fit_proton_resoBothComb[epart], 1, 3,  kGray+4);
    fit_proton_resoBothComb[epart]->Draw("same");

    // DrawGammaSetMarker(h_pion_fhcal_E_allclus[epart], 24, 2, kRed-8, kRed-8);
    // h_pion_fhcal_E_allclus[epart]->Draw("same");
    // DrawGammaSetMarkerTF1(fit_reso_NCell[epart], 1, 2, kRed-4);
    // fit_reso_NCell[epart]->Draw("same");

    DrawGammaSetMarker(h_pion_ecalfhcalcomb_E_calib[epart], 20, 3, kRed+2, kRed+2);
    DrawGammaSetMarkerTF1(fit_resoBothComb_calib[epart], 1, 3,  kRed+4);
    h_pion_ecalfhcalcomb_E_calib[epart]->Draw("same");
    fit_resoBothComb_calib[epart]->Draw("same");

    DrawGammaSetMarker(h_gamma_ecalfhcalcomb_E[epart], 47, 3, kOrange+2, kOrange+2);
    DrawGammaSetMarkerTF1(fit_gamma_resoBothComb[epart], 1, 3,  kOrange+4);
    h_gamma_ecalfhcalcomb_E[epart]->Draw("same");
    fit_gamma_resoBothComb[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV particles",partE[epart]),0.15,0.60,1.0*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("in #eta > 1.2"),0.15,0.55,1.0*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);


    legendHE_proton           = GetAndSetLegend2(0.15, 0.95-(4*textSizeLabelsRel), 0.7, 0.95,1.1*textSizeLabelsPixel, 1, "", 43, 0.15);
    // legendHE_proton->AddEntry(h_pion_fhcal_E_allclus[epart],Form("#pi^{-} : mean: %1.2f, #sigma/#it{E}: %1.2f%%",fit_reso_NCell[epart]->GetParameter(1),100*fit_reso_NCell[epart]->GetParameter(2)/fit_reso_NCell[epart]->GetParameter(1)),"pl");
    // legendHE_proton->AddEntry(h_proton_ecalfhcal_E[epart],Form("HCAL: no fit"),"pl");
    legendHE_proton->AddEntry(h_pion_ecalfhcalcomb_E_calib[epart],Form("#pi^{-}: mean: %1.2f, #sigma/#it{E}: %1.2f%%",fit_resoBothComb_calib[epart]->GetParameter(1),100*fit_resoBothComb_calib[epart]->GetParameter(2)/fit_resoBothComb_calib[epart]->GetParameter(1)),"pl");
    legendHE_proton->AddEntry(h_proton_ecalfhcalcomb_E[epart],Form("proton: mean: %1.2f, #sigma/#it{E}: %1.2f%%",fit_proton_resoBothComb[epart]->GetParameter(1),100*fit_proton_resoBothComb[epart]->GetParameter(2)/fit_proton_resoBothComb[epart]->GetParameter(1)),"pl");
    legendHE_proton->AddEntry(h_gamma_ecalfhcalcomb_E[epart],Form("photon: mean: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothComb[epart]->GetParameter(1),100*fit_gamma_resoBothComb[epart]->GetParameter(2)/fit_gamma_resoBothComb[epart]->GetParameter(1)),"pl");
    // legendNonLin->AddEntry((TObject*)0,"","");
    legendHE_proton->Draw();
    canvasECCE->Print(Form("%s/ECAL_HCAL_resolution_allpart.%s", outputDir.Data(), suffix.Data()));
 }


  split_canvas(cPNGBoth, "cPNGPosition", activeEboth_gamma);
  histEPlotDummyBoth                           = new TH2F("histEPlotDummyBoth","histEPlotDummyBoth",1000,-0.5, 0.5,1000,0, 0.085);
  SetStyleHistoTH2ForGraphs(histEPlotDummyBoth, "#delta#eta, #delta#varphi, #deltaR","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummyBoth->GetXaxis()->SetNoExponent();
  histEPlotDummyBoth->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyBoth->GetXaxis()->SetMoreLogLabels(kTRUE);

  // TF1 *fit_gamma_resoBothComb[nEne] = {NULL};
  // TF1 *fit_gamma_resoBoth[nEne] = {NULL};
  // TF1 *fit_gamma_resoBothEMC[nEne] = {NULL};
  // TLegend* legendHE_gamma = NULL;
  // h_pion__dEta_ecalfhcal_E[epart]->Fill(rec_eta2-rec_eta1);
  //           h_pion__dPhi_ecalfhcal_E[epart]->Fill(rec_phi2-rec_phi1);
  //           h_pion__dR_ecalfhcal_E[epart]->Fill(TMath::Sqrt(TMath::Power(rec_phi2-rec_phi1,2)+TMath::Power(rec_eta2-rec_eta1,2)));
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    if(!h_pion__dR_ecalfhcal_E[epart]) continue;
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    DrawGammaSetMarker(h_pion__dR_ecalfhcal_E[epart], 28, 2, kGreen+2, kGreen+2);
    h_pion__dR_ecalfhcal_E[epart]->Draw("same");
    // fit_gamma_resoBothEMC[epart] = new TF1(Form("fit_gamma_resoBothEMC_%d",epart), "gaus", 0.8,0.9);
    // h_pion__dR_ecalfhcal_E[epart]->Fit(fit_gamma_resoBothEMC[epart],"L0RMEQ");
    // fit_gamma_resoBothEMC[epart]->SetNpx(10000);
    // DrawGammaSetMarkerTF1(fit_gamma_resoBothEMC[epart], 1, 2, kGreen+4);
    // fit_gamma_resoBothEMC[epart]->Draw("same");
    // h_gamma_resolution_fecal_E->SetBinContent(h_gamma_resolution_fecal_E->FindBin(partE[epart]),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1));
    // h_gamma_resolution_fecal_1oE->SetBinContent(h_gamma_resolution_fecal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1));


    DrawGammaSetMarker(h_pion__dPhi_ecalfhcal_E[epart], 27, 2, kBlue+2, kBlue+2);
    h_pion__dPhi_ecalfhcal_E[epart]->Draw("same");
    // fit_gamma_resoBoth[epart] = new TF1(Form("fit_gamma_resoBoth_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);
    // h_pion__dPhi_ecalfhcal_E[epart]->Fit(fit_gamma_resoBoth[epart],"L0RMEQ");
    // fit_gamma_resoBoth[epart]->SetNpx(10000);
    // DrawGammaSetMarkerTF1(fit_gamma_resoBoth[epart], 1, 2, kBlack);
    // fit_gamma_resoBoth[epart]->Draw("same");

    DrawGammaSetMarker(h_pion__dEta_ecalfhcal_E[epart], 20, 2, kRed+2, kRed+2);
    h_pion__dEta_ecalfhcal_E[epart]->Draw("same");
    // fit_gamma_resoBothComb[epart] = new TF1(Form("fit_gamma_resoBothComb_%d",epart), "gaus", 0.8,0.9); //6+0.005*partE[epart]
    // h_pion__dEta_ecalfhcal_E[epart]->Fit(fit_gamma_resoBothComb[epart],"L0RMEQ");
    // fit_gamma_resoBothComb[epart]->SetNpx(10000);
    // DrawGammaSetMarkerTF1(fit_gamma_resoBothComb[epart], 1, 2,  kRed+4);
    // fit_gamma_resoBothComb[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV #pi^{-}",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    // legendHE_gamma           = GetAndSetLegend2(0.20, 0.92-(4*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    // legendHE_gamma->AddEntry(h_gamma_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothEMC[epart]->GetParameter(1),fit_gamma_resoBothEMC[epart]->GetParameter(2),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1)),"pl");
    // legendHE_gamma->AddEntry(h_gamma_ecalfhcal_E[epart],Form("HCAL: no fit"),"pl");
    // // legendHE_gamma->AddEntry(h_gamma_ecalfhcal_E[epart],Form("HCAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBoth[epart]->GetParameter(1),fit_gamma_resoBoth[epart]->GetParameter(0),100*fit_gamma_resoBoth[epart]->GetParameter(0)/fit_gamma_resoBoth[epart]->GetParameter(1)),"pl");
    // legendHE_gamma->AddEntry(h_gamma_ecalfhcalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothComb[epart]->GetParameter(1),fit_gamma_resoBothComb[epart]->GetParameter(2),100*fit_gamma_resoBothComb[epart]->GetParameter(2)/fit_gamma_resoBothComb[epart]->GetParameter(1)),"pl");
    // // legendNonLin->AddEntry((TObject*)0,"","");
    // legendHE_gamma->Draw();
    padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_HCAL_dEtadPhidR_cpion.%s", outputDir.Data(), suffix.Data()));
  padnum =0;
  TF1 *fit_reso_R_HCAL[nEne] = {NULL};
  TF1 *fit_reso_Eta_HCAL[nEne] = {NULL};
  TF1 *fit_reso_Phi_HCAL[nEne] = {NULL};
  TF1 *fit_reso_R_ECAL[nEne] = {NULL};
  TF1 *fit_reso_Eta_ECAL[nEne] = {NULL};
  TF1 *fit_reso_Phi_ECAL[nEne] = {NULL};
  TH1F* h_HCAL_pion_R_reso_E 	    = new TH1F("h_HCAL_pion_R_reso_E", "", 200, 0.0, 100.0);
  TH1F* h_HCAL_pion_Eta_reso_E 	    = new TH1F("h_HCAL_pion_Eta_reso_E", "", 200, 0.0, 100.0);
  TH1F* h_HCAL_pion_Phi_reso_E 	    = new TH1F("h_HCAL_pion_Phi_reso_E", "", 200, 0.0, 100.0);
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    if(!h_pion__dR_hcal_E[epart]) continue;
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    DrawGammaSetMarker(h_pion__dR_hcal_E[epart], 28, 2, kGreen+2, kGreen+2);
    h_pion__dR_hcal_E[epart]->Draw("same");
    fit_reso_R_HCAL[epart] = new TF1(Form("fit_reso_R_HCAL_%d",epart), "gaus", 0.00,0.2);
    h_pion__dR_hcal_E[epart]->Fit(fit_reso_R_HCAL[epart],"L0RMEQ");
    fit_reso_R_HCAL[epart] = new TF1(Form("fit_reso_R_HCAL_%d",epart), "gaus", fit_reso_R_HCAL[epart]->GetParameter(1)-1*fit_reso_R_HCAL[epart]->GetParameter(2),fit_reso_R_HCAL[epart]->GetParameter(1)+1*fit_reso_R_HCAL[epart]->GetParameter(2));
    h_pion__dR_hcal_E[epart]->Fit(fit_reso_R_HCAL[epart],"L0RMEQ");
    fit_reso_R_HCAL[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso_R_HCAL[epart], 1, 2, kGreen+4);
    fit_reso_R_HCAL[epart]->Draw("same");
    // h_HCAL_pion_R_reso_E->SetBinContent(h_HCAL_pion_R_reso_E->FindBin(partE[epart]),100*fit_reso_R_HCAL[epart]->GetParameter(2)/fit_reso_R_HCAL[epart]->GetParameter(1));
    h_HCAL_pion_R_reso_E->SetBinContent(h_HCAL_pion_R_reso_E->FindBin(partE[epart]),fit_reso_R_HCAL[epart]->GetParameter(1));


    DrawGammaSetMarker(h_pion__dPhi_hcal_E[epart], 27, 2, kBlue+2, kBlue+2);
    h_pion__dPhi_hcal_E[epart]->Draw("same");
    fit_reso_Phi_HCAL[epart] = new TF1(Form("fit_reso_Phi_HCAL_%d",epart), "gaus", -0.1,0.2);
    h_pion__dPhi_hcal_E[epart]->Fit(fit_reso_Phi_HCAL[epart],"L0RMEQ");
    h_pion__dPhi_hcal_E[epart]->Draw("same");
    fit_reso_Phi_HCAL[epart] = new TF1(Form("fit_reso_Phi_HCAL_%d",epart), "gaus", fit_reso_Phi_HCAL[epart]->GetParameter(1)-1*fit_reso_Phi_HCAL[epart]->GetParameter(2),fit_reso_Phi_HCAL[epart]->GetParameter(1)+1*fit_reso_Phi_HCAL[epart]->GetParameter(2));
    h_pion__dPhi_hcal_E[epart]->Fit(fit_reso_Phi_HCAL[epart],"L0RMEQ");
    fit_reso_Phi_HCAL[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso_Phi_HCAL[epart], 1, 2, kBlue+4);
    fit_reso_Phi_HCAL[epart]->Draw("same");
    // h_HCAL_pion_Phi_reso_E->SetBinContent(h_HCAL_pion_Phi_reso_E->FindBin(partE[epart]),100*fit_reso_Phi_HCAL[epart]->GetParameter(2)/fit_reso_Phi_HCAL[epart]->GetParameter(1));
    h_HCAL_pion_Phi_reso_E->SetBinContent(h_HCAL_pion_Phi_reso_E->FindBin(partE[epart]),fit_reso_Phi_HCAL[epart]->GetParameter(1));


    DrawGammaSetMarker(h_pion__dEta_hcal_E[epart], 20, 2, kRed+2, kRed+2);
    h_pion__dEta_hcal_E[epart]->Draw("same");
    fit_reso_Eta_HCAL[epart] = new TF1(Form("fit_reso_Eta_HCAL_%d",epart), "gaus", 0.00,0.2);
    h_pion__dEta_hcal_E[epart]->Fit(fit_reso_Eta_HCAL[epart],"L0RMEQ");
    fit_reso_Eta_HCAL[epart] = new TF1(Form("fit_reso_Eta_HCAL_%d",epart), "gaus", fit_reso_Eta_HCAL[epart]->GetParameter(1)-1*fit_reso_Eta_HCAL[epart]->GetParameter(2),fit_reso_Eta_HCAL[epart]->GetParameter(1)+1*fit_reso_Eta_HCAL[epart]->GetParameter(2));
    h_pion__dEta_hcal_E[epart]->Fit(fit_reso_Eta_HCAL[epart],"L0RMEQ");
    fit_reso_Eta_HCAL[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso_Eta_HCAL[epart], 1, 2, kRed+4);
    fit_reso_Eta_HCAL[epart]->Draw("same");
    drawLatexAdd(Form("%1.1f GeV #pi^{-}",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // h_HCAL_pion_Eta_reso_E->SetBinContent(h_HCAL_pion_Eta_reso_E->FindBin(partE[epart]),100*fit_reso_Eta_HCAL[epart]->GetParameter(2)/fit_reso_Eta_HCAL[epart]->GetParameter(1));
    h_HCAL_pion_Eta_reso_E->SetBinContent(h_HCAL_pion_Eta_reso_E->FindBin(partE[epart]),fit_reso_Eta_HCAL[epart]->GetParameter(1));

  legendHE_gamma           = GetAndSetLegend2(0.20, 0.92-(4*textSizeLabelsRel), 0.9, 0.92,0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendHE_gamma->AddEntry(h_pion__dR_hcal_E[epart],"R","p");
  legendHE_gamma->AddEntry(h_pion__dPhi_hcal_E[epart],"#varphi","p");
  legendHE_gamma->AddEntry(h_pion__dEta_hcal_E[epart],"#eta","p");
  legendHE_gamma->Draw();

    // legendHE_gamma           = GetAndSetLegend2(0.20, 0.92-(4*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    // legendHE_gamma->AddEntry(h_gamma_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothEMC[epart]->GetParameter(1),fit_gamma_resoBothEMC[epart]->GetParameter(2),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1)),"pl");
    // legendHE_gamma->AddEntry(h_gamma_hcal_E[epart],Form("HCAL: no fit"),"pl");
    // // legendHE_gamma->AddEntry(h_gamma_hcal_E[epart],Form("HCAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBoth[epart]->GetParameter(1),fit_gamma_resoBoth[epart]->GetParameter(0),100*fit_gamma_resoBoth[epart]->GetParameter(0)/fit_gamma_resoBoth[epart]->GetParameter(1)),"pl");
    // legendHE_gamma->AddEntry(h_gamma_hcalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothComb[epart]->GetParameter(1),fit_gamma_resoBothComb[epart]->GetParameter(2),100*fit_gamma_resoBothComb[epart]->GetParameter(2)/fit_gamma_resoBothComb[epart]->GetParameter(1)),"pl");
    // // legendNonLin->AddEntry((TObject*)0,"","");
    // legendHE_gamma->Draw();
    padnum++;
 }
  cPNGBoth->Print(Form("%s/HCAL_dEtadPhidR_cpion.%s", outputDir.Data(), suffix.Data()));
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    if(!h_pion__dR_ecal_E[epart]) continue;
  	cPNGBoth->cd(padnum+1);
    histEPlotDummyBoth->DrawCopy();
    DrawGammaSetMarker(h_pion__dR_ecal_E[epart], 28, 2, kGreen+2, kGreen+2);
    h_pion__dR_ecal_E[epart]->Draw("same");
    fit_reso_R_ECAL[epart] = new TF1(Form("fit_reso_R_ECAL_%d",epart), "gaus", 0.00,0.2);
    h_pion__dR_ecal_E[epart]->Fit(fit_reso_R_ECAL[epart],"L0RMEQ");
    fit_reso_R_ECAL[epart] = new TF1(Form("fit_reso_R_ECAL_%d",epart), "gaus", fit_reso_R_ECAL[epart]->GetParameter(1)-1*fit_reso_R_ECAL[epart]->GetParameter(2),fit_reso_R_ECAL[epart]->GetParameter(1)+1*fit_reso_R_ECAL[epart]->GetParameter(2));
    h_pion__dR_ecal_E[epart]->Fit(fit_reso_R_ECAL[epart],"L0RMEQ");
    fit_reso_R_ECAL[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso_R_ECAL[epart], 1, 2, kGreen+4);
    fit_reso_R_ECAL[epart]->Draw("same");

    DrawGammaSetMarker(h_pion__dPhi_ecal_E[epart], 27, 2, kBlue+2, kBlue+2);
    h_pion__dPhi_ecal_E[epart]->Draw("same");
    fit_reso_Phi_ECAL[epart] = new TF1(Form("fit_reso_Phi_ECAL_%d",epart), "gaus", 0.0,0.2);
    h_pion__dPhi_ecal_E[epart]->Fit(fit_reso_Phi_ECAL[epart],"L0RMEQ");
    fit_reso_Phi_ECAL[epart] = new TF1(Form("fit_reso_Phi_ECAL_%d",epart), "gaus", fit_reso_Phi_ECAL[epart]->GetParameter(1)-1*fit_reso_Phi_ECAL[epart]->GetParameter(2),fit_reso_Phi_ECAL[epart]->GetParameter(1)+1*fit_reso_Phi_ECAL[epart]->GetParameter(2));
    h_pion__dPhi_ecal_E[epart]->Fit(fit_reso_Phi_ECAL[epart],"L0RMEQ");
    fit_reso_Phi_ECAL[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso_Phi_ECAL[epart], 1, 2, kBlue+4);
    fit_reso_Phi_ECAL[epart]->Draw("same");

    DrawGammaSetMarker(h_pion__dEta_ecal_E[epart], 20, 2, kRed+2, kRed+2);
    h_pion__dEta_ecal_E[epart]->Draw("same");
    fit_reso_Eta_ECAL[epart] = new TF1(Form("fit_reso_Eta_ECAL_%d",epart), "gaus", -0.1,0.1);
    h_pion__dEta_ecal_E[epart]->Fit(fit_reso_Eta_ECAL[epart],"L0RMEQ");
    fit_reso_Eta_ECAL[epart] = new TF1(Form("fit_reso_Eta_ECAL_%d",epart), "gaus", fit_reso_Eta_ECAL[epart]->GetParameter(1)-1*fit_reso_Eta_ECAL[epart]->GetParameter(2),fit_reso_Eta_ECAL[epart]->GetParameter(1)+1*fit_reso_Eta_ECAL[epart]->GetParameter(2));
    h_pion__dEta_ecal_E[epart]->Fit(fit_reso_Eta_ECAL[epart],"L0RMEQ");
    fit_reso_Eta_ECAL[epart]->SetNpx(10000);
    DrawGammaSetMarkerTF1(fit_reso_Eta_ECAL[epart], 1, 2, kRed+4);
    fit_reso_Eta_ECAL[epart]->Draw("same");

    drawLatexAdd(Form("%1.1f GeV #pi^{-}",partE[epart]),0.90,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    // legendHE_gamma           = GetAndSetLegend2(0.20, 0.92-(4*textSizeLabelsRel), 0.9, 0.92,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    // legendHE_gamma->AddEntry(h_gamma_ecal_E[epart],Form("ECAL: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothEMC[epart]->GetParameter(1),fit_gamma_resoBothEMC[epart]->GetParameter(2),100*fit_gamma_resoBothEMC[epart]->GetParameter(2)/fit_gamma_resoBothEMC[epart]->GetParameter(1)),"pl");
    // legendHE_gamma->AddEntry(h_gamma_ecal_E[epart],Form("ecal: no fit"),"pl");
    // // legendHE_gamma->AddEntry(h_gamma_ecal_E[epart],Form("ecal: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBoth[epart]->GetParameter(1),fit_gamma_resoBoth[epart]->GetParameter(0),100*fit_gamma_resoBoth[epart]->GetParameter(0)/fit_gamma_resoBoth[epart]->GetParameter(1)),"pl");
    // legendHE_gamma->AddEntry(h_gamma_ecalcomb_E[epart],Form("Comb: mean: %1.2f, #sigma: %1.2f, #sigma/#it{E}: %1.2f%%",fit_gamma_resoBothComb[epart]->GetParameter(1),fit_gamma_resoBothComb[epart]->GetParameter(2),100*fit_gamma_resoBothComb[epart]->GetParameter(2)/fit_gamma_resoBothComb[epart]->GetParameter(1)),"pl");
    // // legendNonLin->AddEntry((TObject*)0,"","");
    // legendHE_gamma->Draw();
    padnum++;
 }
  cPNGBoth->Print(Form("%s/ECAL_dEtadPhidR_cpion.%s", outputDir.Data(), suffix.Data()));

  TH2F * histEPlotDummyMaterial                           = new TH2F("histEPlotDummyMaterial","histEPlotDummyMaterial",1000,0, 1.5,1000,0, 0.09);
  SetStyleHistoTH2ForGraphs(histEPlotDummyMaterial, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  histEPlotDummyMaterial->GetXaxis()->SetNoExponent();
  histEPlotDummyMaterial->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyMaterial->GetXaxis()->SetMoreLogLabels(kTRUE);

  split_canvas(cPNG, "cPNG4", activeEboth);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart] && !useEBothCPion[epart]) continue;
    if(!useEBothCPion[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyMaterial->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV #pi^{-}}",partE[epart]),0.90,0.26,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart] && h_pion_fhcal_E[epart]){
      DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      h_pion_fhcal_E[epart]->Draw("same");
      DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      fit_reso[epart]->Draw("same");
      drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    if(useEBothCPion[epart] && h_pion_ecalfhcal_E[epart]){
      DrawGammaSetMarker(h_pion_ecalfhcal_E[epart], 20, 2, kAzure+6, kAzure+6);
      h_pion_ecalfhcal_E[epart]->Draw("same");

      DrawGammaSetMarkerTF1(fit_resoBoth[epart], 1, 2, kBlack);
      fit_resoBoth[epart]->Draw("same");

      drawLatexAdd(Form("mean: %1.2f ",fit_resoBoth[epart]->GetParameter(1)),0.9,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("sigma: %1.2f ",fit_resoBoth[epart]->GetParameter(2)),0.9,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1)),0.9,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#color[866]{#font[62]{FHCAL + FECAL}}",0.9,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    }
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_pion_material_resolution.%s", outputDir.Data(), suffix.Data()));

  TH2F * histEPlotDummyEHEE                           = new TH2F("histEPlotDummyEHEE","histEPlotDummyEHEE",2000,0.0, 200,2000,0.0, 200);
  SetStyleHistoTH2ForGraphs(histEPlotDummyEHEE, "#it{E}_{ECAL}","#it{E}_{HCAL}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  histEPlotDummyEHEE->GetXaxis()->SetNoExponent();
  histEPlotDummyEHEE->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyEHEE->GetXaxis()->SetMoreLogLabels(kTRUE);

  split_canvas(cPNG, "cPNG4", activeEboth);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyEHEE->GetXaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyEHEE->GetYaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyEHEE->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV #pi^{-}}",partE[epart]),0.90,0.86,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useE[epart]){
      // DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      // h2D_pion_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0, 100, 0.0, 2.0);
      h2D_pion_ecalfhcal_E[epart]->Draw("colz,same");
      // DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }
    // if(useEBothCPion[epart]){
    //   DrawGammaSetMarker(h_pion_ecalfhcal_E[epart], 20, 2, kAzure+6, kAzure+6);
    //   h_pion_ecalfhcal_E[epart]->Draw("same");

    //   DrawGammaSetMarkerTF1(fit_resoBoth[epart], 1, 2, kBlack);
    //   fit_resoBoth[epart]->Draw("same");

    //   drawLatexAdd(Form("mean: %1.2f ",fit_resoBoth[epart]->GetParameter(1)),0.9,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    //   drawLatexAdd(Form("sigma: %1.2f ",fit_resoBoth[epart]->GetParameter(2)),0.9,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    //   drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_resoBoth[epart]->GetParameter(2)/fit_resoBoth[epart]->GetParameter(1)),0.9,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    //   drawLatexAdd("#color[866]{#font[62]{FHCAL + FECAL}}",0.9,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // }
    padnum++;
 }
  cPNG->Print(Form("%s/Energy_HCAL_ECAL_cpion.%s", outputDir.Data(), suffix.Data()));

  split_canvas(cPNG, "cPNG4", activeEboth_gamma);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    if(!h2D_gamma_ecalfhcal_E[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyEHEE->GetXaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyEHEE->GetYaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyEHEE->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV #gamma}",partE[epart]),0.90,0.86,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useEBothgamma[epart]){
      // DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      // h2D_pion_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0, 100, 0.0, 2.0);
      h2D_gamma_ecalfhcal_E[epart]->Draw("colz,same");
      // DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }

    padnum++;
 }
  cPNG->Print(Form("%s/Energy_HCAL_ECAL_gamma.%s", outputDir.Data(), suffix.Data()));
  split_canvas(cPNG, "cPNG4", activeEboth_gamma);

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothproton[epart]) continue;
    if(!h2D_proton_ecalfhcal_E[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyEHEE->GetXaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyEHEE->GetYaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyEHEE->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV proton}",partE[epart]),0.90,0.86,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useEBothproton[epart]){
      // DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      // h2D_pion_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0, 100, 0.0, 2.0);
      h2D_proton_ecalfhcal_E[epart]->Draw("colz,same");
      // DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }

    padnum++;
 }
  cPNG->Print(Form("%s/Energy_HCAL_ECAL_proton.%s", outputDir.Data(), suffix.Data()));


  TH2F * histEPlotDummyENCells                           = new TH2F("histEPlotDummyENCells","histEPlotDummyENCells",2000,0.0, 200,2000,0.0, 200);
  SetStyleHistoTH2ForGraphs(histEPlotDummyENCells, "#it{E}_{HCAL}^{cluster}","#it{N}_{towers}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.9);
  histEPlotDummyENCells->GetXaxis()->SetNoExponent();
  histEPlotDummyENCells->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummyENCells->GetXaxis()->SetMoreLogLabels(kTRUE);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothgamma[epart]) continue;
    if(!h2D_NTowers_HCAL_gamma[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyENCells->GetXaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyENCells->GetYaxis()->SetRangeUser(0,30);
    histEPlotDummyENCells->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV #pi^{-}}",partE[epart]),0.90,0.86,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useEBothgamma[epart]){
      // DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      // h2D_pion_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0, 100, 0.0, 2.0);
      h2D_NTowers_HCAL_gamma[epart]->Draw("colz,same");
      // DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }

    padnum++;
 }
  cPNG->Print(Form("%s/NCells_HCAL_ECAL_gamma.%s", outputDir.Data(), suffix.Data()));

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothCPion[epart]) continue;
    if(!h2D_NTowers_HCAL_pi[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyENCells->GetXaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyENCells->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV #pi^{-}}",partE[epart]),0.90,0.86,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useEBothCPion[epart]){
      // DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      // h2D_pion_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0, 100, 0.0, 2.0);
      h2D_NTowers_HCAL_pi[epart]->Draw("colz,same");
      // DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }

    padnum++;
 }
  cPNG->Print(Form("%s/NCells_HCAL_ECAL_pion.%s", outputDir.Data(), suffix.Data()));

  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useEBothproton[epart]) continue;
    if(!h2D_NTowers_HCAL_proton[epart]) continue;
    // if(partE[epart]==70 || partE[epart]==80 ) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummyENCells->GetXaxis()->SetRangeUser(0,partE[epart]*1.15);
    histEPlotDummyENCells->DrawCopy();
    drawLatexAdd(Form("#font[62]{%1.1f GeV proton}",partE[epart]),0.90,0.86,1.5*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("rec. in FHCAL",0.90,0.32,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    if(useEBothproton[epart]){
      // DrawGammaSetMarker(h_pion_fhcal_E[epart], 20, 2, kBlue+2, kBlue+2);
      // h2D_pion_ecalfhcal_E[epart] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0, 100, 0.0, 2.0);
      h2D_NTowers_HCAL_proton[epart]->Draw("colz,same");
      // DrawGammaSetMarkerTF1(fit_reso[epart], 1, 2, kBlack);
      // fit_reso[epart]->Draw("same");
      // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.20,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.20,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.20,0.73,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd("#color[602]{#font[62]{FHCAL}}",0.25,0.65,1.5*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    }

    padnum++;
 }
  cPNG->Print(Form("%s/NCells_HCAL_ECAL_proton.%s", outputDir.Data(), suffix.Data()));


// RESOLUTION PLOT

  TCanvas* cReso = new TCanvas("cReso","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso, 0.095, 0.01, 0.01, 0.105);
  TH2F * histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 89,1000,0, 39);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
    DrawGammaSetMarker(h_pion_resolution_fhcal_E, 46, 3, kBlack, kBlack);

  h_pion_resolution_fhcal_E->Draw("same,p");
  TF1* fit_reso_E = new TF1("fit_reso_E", "[0]+([1]/TMath::Sqrt(x))+([2]/x)", 0,200);
  fit_reso_E->SetParLimits(0,0.,15);
  // fit_reso_E->SetParLimits(0,9.,10.);
  // fit_reso_E->SetParLimits(1,30.,90);
  fit_reso_E->FixParameter(2,0.08);
  h_pion_resolution_fhcal_E->Fit(fit_reso_E,"QRMNE");
  fit_reso_E->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_E, 3, 3, kRed+2);
  fit_reso_E->Draw("same");
  drawLatexAdd("#pi^{-} in FHCAL",0.90,0.85,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_E->GetParameter(1),fit_reso_E->GetParameter(0)),0.90,0.75,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/HCAL_resolution_Fitted.%s", outputDir.Data(), suffix.Data()));

  histResoDummy->Draw();

  DrawGammaSetMarker(h_pion_resolution_fhcal_E, 46, 3, kBlue+2, kBlue+2);
  h_pion_resolution_fhcal_E->SetLineStyle(3);
  h_pion_resolution_fhcal_E->SetLineWidth(3);
  h_pion_resolution_fhcal_E->Draw("same,p");
  DrawGammaSetMarkerTF1(fit_reso_E, 3, 3, kBlue+2);
  fit_reso_E->Draw("same");

  DrawGammaSetMarker(h_pion_resolution_fecalfhcal_E, 24, 3, kBlack, kBlack);
  h_pion_resolution_fecalfhcal_E->SetLineStyle(5);
  h_pion_resolution_fecalfhcal_E->SetLineWidth(3);
  h_pion_resolution_fecalfhcal_E->Draw("same,p");
  TF1* fit_resoboth_E = new TF1("fit_resoboth_E", "[0]+[1]/TMath::Sqrt(x)", 0,200);
  fit_resoboth_E->SetParLimits(0,0.,15);
  h_pion_resolution_fecalfhcal_E->Fit(fit_resoboth_E,"QRMNE+");
  fit_resoboth_E->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_resoboth_E, 5, 3, kBlack);
  fit_resoboth_E->Draw("same");

  DrawGammaSetMarker(h_gamma_resolution_fhcal_E, 27, 3, kGreen+2, kGreen+2);
  h_gamma_resolution_fhcal_E->SetLineStyle(4);
  h_gamma_resolution_fhcal_E->SetLineWidth(3);
  h_gamma_resolution_fhcal_E->Draw("same,p");
  TF1* fit_resogamma_E = new TF1("fit_resogamma_E", "[0]+[1]/TMath::Sqrt(x)", 0,200);
  fit_resogamma_E->SetParLimits(0,0.,15);
  h_gamma_resolution_fhcal_E->Fit(fit_resogamma_E,"QRMNE+");
  fit_resogamma_E->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_resogamma_E, 4, 3, kGreen+2);
  fit_resogamma_E->Draw("same");

  TLegend* legendNonLin           = GetAndSetLegend2(0.25, 0.95-(5*textSizeLabelsRel), 0.9, 0.95,textSizeLabelsPixel, 1, "", 43, 0.15);
  legendNonLin->AddEntry(h_pion_resolution_fhcal_E,Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_E->GetParameter(1),fit_reso_E->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_pion_resolution_fecalfhcal_E,Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_E->GetParameter(1),fit_resoboth_E->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_gamma_resolution_fhcal_E,Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resogamma_E->GetParameter(1),fit_resogamma_E->GetParameter(0)),"pl");
  // legendNonLin->AddEntry((TObject*)0,"","");
  legendNonLin->Draw();

  // drawLatexAdd(Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_E->GetParameter(1),fit_reso_E->GetParameter(0)),0.93,0.90,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // drawLatexAdd(Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_E->GetParameter(1),fit_resoboth_E->GetParameter(0)),0.93,0.80,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  // drawLatexAdd(Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resogamma_E->GetParameter(1),fit_resogamma_E->GetParameter(0)),0.93,0.70,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  cReso->Print(Form("%s/FECALFHCAL_resolution_Fitted.%s", outputDir.Data(), suffix.Data()));

// RESOLUTION PLOT 2

  TCanvas* cReso2x = new TCanvas("cReso2x","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2x, 0.1, 0.01, 0.01, 0.105);
  cReso2x->SetLogx();
  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,1, 159,1000,0.31, 1.05);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","< #it{E}_{rec} / #it{E}_{true} >", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
  // DrawGammaSetMarker(h_gamma_hcal_mean_E, 20, 3, kOrange-8, kOrange-8);
  DrawGammaSetMarker(h_gamma_ecal_mean_E, 20, 3, kOrange-8, kOrange-8);
  DrawGammaSetMarker(h_gamma_ecal_mean_E_calibrated, 28, 3, kOrange+2, kOrange+2);
  DrawGammaSetMarker(h_cpion_comb_mean_E, 30, 3, kRed-8, kRed-8);
  DrawGammaSetMarker(h_cpion_comb_mean_E_calib, 24, 3, kRed+2, kRed+2);
  DrawGammaSetMarker(h_cpion_comb_mean_E_calibrated, 46, 3, kRed, kRed);
  DrawGammaSetMarker(h_proton_comb_mean_E, 28, 3, kMagenta-8, kMagenta-8);
  DrawGammaSetMarker(h_proton_comb_mean_E_calibrated, 27, 3, kMagenta+2, kMagenta+2);

  // h_gamma_hcal_mean_E->Draw("same,p");
  h_gamma_ecal_mean_E->Draw("same,p");
  h_gamma_ecal_mean_E_calibrated->Draw("same,p");
  h_cpion_comb_mean_E->Draw("same,p");
  h_cpion_comb_mean_E_calib->Draw("same,p");
  h_cpion_comb_mean_E_calibrated->Draw("same,p");
  h_proton_comb_mean_E->Draw("same,p");
  h_proton_comb_mean_E_calibrated->Draw("same,p");

  TF1* Function_Calib_ECAL   = new TF1("Function_Calib_ECAL",   "( [0]  )",   2, 220);
  h_gamma_ecal_mean_E->Fit(Function_Calib_ECAL,"NRME+","",1.0, 250);
  Function_Calib_ECAL->Draw("same");

  // h_gamma_ecal_mean_E->SetBinContent(h_gamma_ecal_mean_E->FindBin(190),0.84);
  // Double_t paramsTB150MeV[5]= {2.65428e+00,4.17288e-01,1.91848e+00,-5.35529e+02,6.82686e+02}; //noshift
  Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
  TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   2, 220);
  FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
  h_cpion_comb_mean_E_calib->Fit(FunctionNL_PowExp150MeV,"NRME+","",1.0, 250);
  FunctionNL_PowExp150MeV->Draw("same");

  TLegend* legendEmean           = GetAndSetLegend2(0.6, 0.15, 0.95, 0.15+(6*textSizeLabelsRel),0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  // legendEmean->AddEntry(h_gamma_hcal_mean_E,"photon (HCAL)","p");
  legendEmean->AddEntry(h_gamma_ecal_mean_E,"photon (ECAL)","p");
  legendEmean->AddEntry(h_gamma_ecal_mean_E_calibrated,"photon (ECAL calib.)","p");
  legendEmean->AddEntry(h_cpion_comb_mean_E,"charged pion","p");
  legendEmean->AddEntry(h_cpion_comb_mean_E_calib,"charged pion (ECAL calib.)","p");
  legendEmean->AddEntry(h_cpion_comb_mean_E_calibrated,"charged pion (full calib.)","p");
  legendEmean->AddEntry(h_proton_comb_mean_E,"proton","p");
  legendEmean->AddEntry(h_proton_comb_mean_E_calibrated,"proton (full calib.)","p");
  legendEmean->Draw();
  cReso2x->Print(Form("%s/MeanE_HCAL_ECAL.%s", outputDir.Data(), suffix.Data()));

  histResoDummy->Draw();
  // DrawGammaSetMarker(h_gamma_hcal_mean_E, 20, 3, kOrange-8, kOrange-8);
  DrawGammaSetMarker(h_gamma_ecal_mean_E, 20, 3, kOrange-8, kOrange-8);
  DrawGammaSetMarker(h_gamma_ecal_mean_E_calibrated, 28, 3, kOrange+2, kOrange+2);
  DrawGammaSetMarker(h_cpion_comb_mean_E, 30, 3, kRed-8, kRed-8);
  DrawGammaSetMarker(h_cpion_comb_mean_E_calib, 24, 3, kRed+2, kRed+2);
  DrawGammaSetMarker(h_cpion_comb_mean_E_calibrated, 46, 3, kRed, kRed);

  // h_gamma_hcal_mean_E->Draw("same,p");
  h_gamma_ecal_mean_E->Draw("same,p");
  h_gamma_ecal_mean_E_calibrated->Draw("same,p");
  h_cpion_comb_mean_E->Draw("same,p");
  h_cpion_comb_mean_E_calib->Draw("same,p");
  h_cpion_comb_mean_E_calibrated->Draw("same,p");

  Function_Calib_ECAL->Draw("same");

  FunctionNL_PowExp150MeV->Draw("same");

  legendEmean           = GetAndSetLegend2(0.6, 0.15, 0.95, 0.15+(5*textSizeLabelsRel),0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  // legendEmean->AddEntry(h_gamma_hcal_mean_E,"photon (HCAL)","p");
  legendEmean->AddEntry(h_gamma_ecal_mean_E,"photon (ECAL)","p");
  legendEmean->AddEntry(h_gamma_ecal_mean_E_calibrated,"photon (ECAL calib.)","p");
  legendEmean->AddEntry(h_cpion_comb_mean_E,"charged pion","p");
  legendEmean->AddEntry(h_cpion_comb_mean_E_calib,"charged pion (ECAL calib.)","p");
  legendEmean->AddEntry(h_cpion_comb_mean_E_calibrated,"charged pion (full calib.)","p");
  legendEmean->Draw();
  cReso2x->Print(Form("%s/MeanE_HCAL_ECAL_calib.%s", outputDir.Data(), suffix.Data()));

  cReso2x = new TCanvas("cReso2x","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2x, 0.12, 0.01, 0.01, 0.105);
  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 69,1000,0.0, 0.19);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","#sigma (rec - true)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.15);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
  DrawGammaSetMarker(h_HCAL_pion_R_reso_E, 28, 3, kGreen+2, kGreen+2);
  DrawGammaSetMarker(h_HCAL_pion_Phi_reso_E, 30, 3, kBlue+2, kBlue+2);
  DrawGammaSetMarker(h_HCAL_pion_Eta_reso_E, 28, 3, kRed+2, kRed+2);

  h_HCAL_pion_R_reso_E->Draw("same,p");
  h_HCAL_pion_Phi_reso_E->Draw("same,p");
  h_HCAL_pion_Eta_reso_E->Draw("same,p");
  TLegend* legendPosReso           = GetAndSetLegend2(0.6, 0.9-(3*textSizeLabelsRel), 0.95, 0.9,0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendPosReso->AddEntry(h_HCAL_pion_R_reso_E,"R","p");
  legendPosReso->AddEntry(h_HCAL_pion_Phi_reso_E,"#varphi","p");
  legendPosReso->AddEntry(h_HCAL_pion_Eta_reso_E,"#eta","p");
  legendPosReso->Draw();
  drawLatexAdd("#pi^{-} in FHCAL",0.90,0.55,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso2x->Print(Form("%s/PositionReso_HCAL.%s", outputDir.Data(), suffix.Data()));

// RESOLUTION PLOT 2

  TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 34);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
    DrawGammaSetMarker(h_pion_resolution_fhcal_1oE, 46, 3, kBlack, kBlack);

  h_pion_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_reso_1oE = new TF1("fit_reso_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_1oE->SetParameter(2,0.3);
  // fit_reso_1oE->SetParLimits(2,0.0,10);
  fit_reso_1oE->FixParameter(2,0.08);
  h_pion_resolution_fhcal_1oE->Fit(fit_reso_1oE,"QRMNE");
  fit_reso_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_1oE, 3, 3, kRed+2);
  fit_reso_1oE->Draw("same");
  drawLatexAdd("#pi^{-} in FHCAL",0.90,0.25,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd(Form("#frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),0.90,0.35,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


  cReso2->Print(Form("%s/HCAL_resolution_1oE_Fitted.%s", outputDir.Data(), suffix.Data()));

  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 45);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
  DrawGammaSetMarker(h_pion_resolution_fhcal_1oE, 46, 3, kBlue+2, kBlue+2);
  h_pion_resolution_fhcal_1oE->SetLineStyle(3);
  h_pion_resolution_fhcal_1oE->SetLineWidth(3);
  h_pion_resolution_fhcal_1oE->Draw("same,p");
  fit_reso_1oE = new TF1("fit_reso_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_1oE->SetParLimits(0,0.,15.);
  fit_reso_1oE->FixParameter(2,0.08);
  h_pion_resolution_fhcal_1oE->Fit(fit_reso_1oE,"QRMNE");
  fit_reso_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_1oE, 3, 3, kBlue+2);
  fit_reso_1oE->Draw("same");

  DrawGammaSetMarker(h_pion_resolution_fecalfhcal_1oE, 24, 3, kBlack, kBlack);
  h_pion_resolution_fhcal_1oE->SetLineStyle(5);
  h_pion_resolution_fhcal_1oE->SetLineWidth(3);
  h_pion_resolution_fecalfhcal_1oE->Draw("same,p");
  TF1* fit_resohcalmat_1oE = new TF1("fit_resohcalmat_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_resohcalmat_1oE->SetParLimits(0,0.,15.);
  fit_resohcalmat_1oE->FixParameter(2,0.08);
  h_pion_resolution_fecalfhcal_1oE->Fit(fit_resohcalmat_1oE,"QRMNE");
  fit_resohcalmat_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_resohcalmat_1oE, 3, 3, kBlack);
  fit_resohcalmat_1oE->Draw("same");

  DrawGammaSetMarker(h_pion_resolution_comb_1oE, 30, 3, kRed+2, kRed+2);
  h_pion_resolution_comb_1oE->SetLineStyle(5);
  h_pion_resolution_comb_1oE->SetLineWidth(3);
  h_pion_resolution_comb_1oE->Draw("same,p");
  TF1* fit_resoboth_1oE = new TF1("fit_resoboth_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_resoboth_1oE->SetParLimits(0,0.,15.);
  fit_resoboth_1oE->FixParameter(2,0.08);
  h_pion_resolution_comb_1oE->Fit(fit_resoboth_1oE,"QRMNE");
  fit_resoboth_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_resoboth_1oE, 3, 3, kRed+2);
  fit_resoboth_1oE->Draw("same");

  DrawGammaSetMarker(h_gamma_resolution_fhcal_1oE, 27, 3, kGreen+2, kGreen+2);
  h_gamma_resolution_fhcal_1oE->SetLineStyle(7);
  h_gamma_resolution_fhcal_1oE->SetLineWidth(3);
  h_gamma_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_reso_gamma_1oE = new TF1("fit_reso_gamma_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_gamma_1oE->SetParLimits(0,0.,15.);
  fit_reso_gamma_1oE->FixParameter(2,0.08);
  h_gamma_resolution_fhcal_1oE->Fit(fit_reso_gamma_1oE,"QRMNE");
  fit_reso_gamma_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_gamma_1oE, 3, 3, kGreen+2);
  fit_reso_gamma_1oE->Draw("same");

  DrawGammaSetMarker(h_gamma_resolution_fecal_1oE, 28, 3, kOrange+2, kOrange+2);
  h_gamma_resolution_fecal_1oE->SetLineStyle(7);
  h_gamma_resolution_fecal_1oE->SetLineWidth(3);
  h_gamma_resolution_fecal_1oE->Draw("same,p");
  TF1* fit_reso_gamma_EMC_1oE = new TF1("fit_reso_gamma_EMC_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_gamma_EMC_1oE->SetParLimits(0,0.,15.);
  fit_reso_gamma_EMC_1oE->FixParameter(2,0.08);
  h_gamma_resolution_fecal_1oE->Fit(fit_reso_gamma_EMC_1oE,"QRMNE");
  fit_reso_gamma_EMC_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_gamma_EMC_1oE, 3, 3, kOrange+2);
  fit_reso_gamma_EMC_1oE->Draw("same");

  DrawGammaSetMarker(h_proton_resolution_fecalfhcal_1oE, 28, 3, kMagenta+2, kMagenta+2);
  h_proton_resolution_fecalfhcal_1oE->SetLineStyle(7);
  h_proton_resolution_fecalfhcal_1oE->SetLineWidth(3);
  h_proton_resolution_fecalfhcal_1oE->Draw("same,p");
  TF1* fit_reso_proton_both_1oE = new TF1("fit_reso_proton_both_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_reso_proton_both_1oE->SetParLimits(0,0.,15.);
  fit_reso_proton_both_1oE->FixParameter(2,0.08);
  h_proton_resolution_fecalfhcal_1oE->Fit(fit_reso_proton_both_1oE,"QRMNE");
  fit_reso_proton_both_1oE->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_reso_proton_both_1oE, 3, 3, kMagenta+2);
  fit_reso_proton_both_1oE->Draw("same");

  // drawLatexAdd(Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resohcalmat_1oE->GetParameter(1),fit_resohcalmat_1oE->GetParameter(0)),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_1oE->GetParameter(1),fit_reso_gamma_1oE->GetParameter(0)),0.15,0.70,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  legendNonLin           = GetAndSetLegend2(0.15, 0.95-(9*textSizeLabelsRel), 0.7, 0.95,0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendNonLin->AddEntry(h_pion_resolution_fhcal_1oE,Form("#pi^{-} in HCAL (w/o mat): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_pion_resolution_fecalfhcal_1oE,Form("#pi^{-} in HCAL (w/ mat): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resohcalmat_1oE->GetParameter(1),fit_resohcalmat_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_pion_resolution_comb_1oE,Form("#pi^{-} in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_1oE->GetParameter(1),fit_resoboth_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_proton_resolution_fecalfhcal_1oE,Form("p in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_proton_both_1oE->GetParameter(1),fit_reso_proton_both_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_gamma_resolution_fhcal_1oE,Form("#gamma in HCAL (w/o mat): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_1oE->GetParameter(1),fit_reso_gamma_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_gamma_resolution_fecal_1oE,Form("#gamma in ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_EMC_1oE->GetParameter(1),fit_reso_gamma_EMC_1oE->GetParameter(0)),"pl");
  // legendNonLin->AddEntry((TObject*)0,"","");
  legendNonLin->Draw();
  cReso2->Print(Form("%s/HCAL_resolution_all_1oE_Fitted.%s", outputDir.Data(), suffix.Data()));

  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 45);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
  DrawGammaSetMarker(h_proton_resolution_fecalfhcal_1oE, 46, 3, kGreen+2, kGreen+2);
  h_proton_resolution_fecalfhcal_1oE->Draw("same,p");
  DrawGammaSetMarkerTF1(fit_reso_proton_both_1oE, 3, 3, kGreen+2);
  fit_reso_proton_both_1oE->Draw("same");

  DrawGammaSetMarker(h_pion_resolution_fhcal_1oE, 24, 3, kRed-8, kRed-8);
  h_pion_resolution_fhcal_1oE->Draw("same,p");
  DrawGammaSetMarkerTF1(fit_reso_1oE, 3, 3, kRed-8);
  fit_reso_1oE->Draw("same");

  DrawGammaSetMarker(h_pion_resolution_comb_1oE, 20, 3, kRed+2, kRed+2);
  h_pion_resolution_comb_1oE->Draw("same,p");
  DrawGammaSetMarkerTF1(fit_resoboth_1oE, 3, 3, kRed+2);
  fit_resoboth_1oE->Draw("same");

  DrawGammaSetMarker(h_gamma_resolution_fecal_1oE, 33, 3, kOrange+2, kOrange+2);
  h_gamma_resolution_fecal_1oE->Draw("same,p");
  DrawGammaSetMarkerTF1(fit_reso_gamma_EMC_1oE, 3, 3, kOrange+2);
  fit_reso_gamma_EMC_1oE->Draw("same");

  // drawLatexAdd(Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resohcalmat_1oE->GetParameter(1),fit_resohcalmat_1oE->GetParameter(0)),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_1oE->GetParameter(1),fit_reso_gamma_1oE->GetParameter(0)),0.15,0.70,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  legendNonLin           = GetAndSetLegend2(0.15, 0.95-(6*textSizeLabelsRel), 0.7, 0.95,0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendNonLin->AddEntry(h_pion_resolution_fhcal_1oE,Form("#pi^{-} in HCAL (w/o mat): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),"pl");
  // legendNonLin->AddEntry(h_pion_resolution_fecalfhcal_1oE,Form("#pi^{-} in HCAL (w/ mat): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resohcalmat_1oE->GetParameter(1),fit_resohcalmat_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_pion_resolution_comb_1oE,Form("#pi^{-} in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_1oE->GetParameter(1),fit_resoboth_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_proton_resolution_fecalfhcal_1oE,Form("p in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_proton_both_1oE->GetParameter(1),fit_reso_proton_both_1oE->GetParameter(0)),"pl");
  // legendNonLin->AddEntry(h_gamma_resolution_fhcal_1oE,Form("#gamma in HCAL (w/o mat): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_1oE->GetParameter(1),fit_reso_gamma_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_gamma_resolution_fecal_1oE,Form("#gamma in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_EMC_1oE->GetParameter(1),fit_reso_gamma_EMC_1oE->GetParameter(0)),"pl");
  // legendNonLin->AddEntry((TObject*)0,"","");
  legendNonLin->Draw();
  cReso2->Print(Form("%s/HCAL_resolution_ECCE_1oE_Fitted.%s", outputDir.Data(), suffix.Data()));

  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 45);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();

  h_pion_resolution_comb_1oE->Draw("same,p");
  fit_resoboth_1oE->Draw("same");

  h_gamma_resolution_fecal_1oE->Draw("same,p");
  fit_reso_gamma_EMC_1oE->Draw("same");

  h_proton_resolution_fecalfhcal_1oE->Draw("same,p");
  fit_reso_proton_both_1oE->Draw("same");

  legendNonLin           = GetAndSetLegend2(0.15, 0.95-(5*textSizeLabelsRel), 0.7, 0.95,0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendNonLin->AddEntry(h_pion_resolution_comb_1oE,Form("#pi^{-} in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_resoboth_1oE->GetParameter(1),fit_resoboth_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_proton_resolution_fecalfhcal_1oE,Form("p in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_proton_both_1oE->GetParameter(1),fit_reso_proton_both_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_gamma_resolution_fecal_1oE,Form("#gamma in HCAL+ECAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_gamma_EMC_1oE->GetParameter(1),fit_reso_gamma_EMC_1oE->GetParameter(0)),"pl");
  legendNonLin->Draw();
  cReso2->Print(Form("%s/HCALECAL_resolution_all_1oE_Fitted.%s", outputDir.Data(), suffix.Data()));

  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 55);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  histResoDummy->Draw();
  DrawGammaSetMarker(h_pion_resolution_fhcal_1oE, 46, 3, kBlue+2, kBlue+2);
  h_pion_resolution_fhcal_1oE->SetLineStyle(3);
  h_pion_resolution_fhcal_1oE->SetLineWidth(3);
  h_pion_resolution_fhcal_1oE->Draw("same,p");
  DrawGammaSetMarkerTF1(fit_reso_1oE, 3, 3, kBlue+2);
  fit_reso_1oE->Draw("same");

  DrawGammaSetMarker(h_cpion_eta1_resolution_fhcal_1oE, 24, 3, kBlack, kBlack);
  h_cpion_eta1_resolution_fhcal_1oE->SetLineStyle(5);
  h_cpion_eta1_resolution_fhcal_1oE->SetLineWidth(3);
  h_cpion_eta1_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_cpion_eta1_hcal = new TF1("fit_cpion_eta1_hcal", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_cpion_eta1_hcal->SetParLimits(0,0.,15.);
  fit_cpion_eta1_hcal->FixParameter(2,0.08);
  h_cpion_eta1_resolution_fhcal_1oE->Fit(fit_cpion_eta1_hcal,"QRMNE");
  fit_cpion_eta1_hcal->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_cpion_eta1_hcal, 3, 3, kBlack);
  fit_cpion_eta1_hcal->Draw("same");

  DrawGammaSetMarker(h_cpion_eta2_resolution_fhcal_1oE, 30, 3, kRed+2, kRed+2);
  h_cpion_eta2_resolution_fhcal_1oE->SetLineStyle(5);
  h_cpion_eta2_resolution_fhcal_1oE->SetLineWidth(3);
  h_cpion_eta2_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_cpion_eta2_hcal = new TF1("fit_cpion_eta2_hcal", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_cpion_eta2_hcal->SetParLimits(0,0.,15.);
  fit_cpion_eta2_hcal->FixParameter(2,0.08);
  h_cpion_eta2_resolution_fhcal_1oE->Fit(fit_cpion_eta2_hcal,"QRMNE");
  fit_cpion_eta2_hcal->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_cpion_eta2_hcal, 3, 3, kRed+2);
  fit_cpion_eta2_hcal->Draw("same");

  DrawGammaSetMarker(h_cpion_eta3_resolution_fhcal_1oE, 27, 3, kGreen+2, kGreen+2);
  h_cpion_eta3_resolution_fhcal_1oE->SetLineStyle(7);
  h_cpion_eta3_resolution_fhcal_1oE->SetLineWidth(3);
  h_cpion_eta3_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_cpion_eta3_hcal = new TF1("fit_cpion_eta3_hcal", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_cpion_eta3_hcal->SetParLimits(0,0.,15.);
  fit_cpion_eta3_hcal->FixParameter(2,0.08);
  h_cpion_eta3_resolution_fhcal_1oE->Fit(fit_cpion_eta3_hcal,"QRMNE");
  fit_cpion_eta3_hcal->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_cpion_eta3_hcal, 3, 3, kGreen+2);
  fit_cpion_eta3_hcal->Draw("same");

  DrawGammaSetMarker(h_cpion_eta4_resolution_fhcal_1oE, 28, 3, kOrange+2, kOrange+2);
  h_cpion_eta4_resolution_fhcal_1oE->SetLineStyle(7);
  h_cpion_eta4_resolution_fhcal_1oE->SetLineWidth(3);
  h_cpion_eta4_resolution_fhcal_1oE->Draw("same,p");
  TF1* fit_cpion_eta4_hcal = new TF1("fit_cpion_eta4_hcal", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
  fit_cpion_eta4_hcal->SetParLimits(0,0.,15.);
  fit_cpion_eta4_hcal->FixParameter(2,0.08);
  h_cpion_eta4_resolution_fhcal_1oE->Fit(fit_cpion_eta4_hcal,"QRMNE");
  fit_cpion_eta4_hcal->SetNpx(10000);
  DrawGammaSetMarkerTF1(fit_cpion_eta4_hcal, 3, 3, kOrange+2);
  fit_cpion_eta4_hcal->Draw("same");

  // drawLatexAdd(Form("#pi^{-} in FHCAL(+FEMCAL): #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_cpion_eta1_hcal->GetParameter(1),fit_cpion_eta1_hcal->GetParameter(0)),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#pi^{-} in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),0.15,0.80,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  // drawLatexAdd(Form("#gamma in FHCAL: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_cpion_eta3_hcal->GetParameter(1),fit_cpion_eta3_hcal->GetParameter(0)),0.15,0.70,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

  legendNonLin           = GetAndSetLegend2(0.15, 0.95-(8*textSizeLabelsRel), 0.7, 0.95,0.9*textSizeLabelsPixel, 1, "", 43, 0.15);
  legendNonLin->AddEntry(h_pion_resolution_fhcal_1oE,Form("1.4<#eta<3.5: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_1oE->GetParameter(1),fit_reso_1oE->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_cpion_eta1_resolution_fhcal_1oE,Form("1.2<#eta<1.7: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_cpion_eta1_hcal->GetParameter(1),fit_cpion_eta1_hcal->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_cpion_eta2_resolution_fhcal_1oE,Form("1.7<#eta<2.5: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_cpion_eta2_hcal->GetParameter(1),fit_cpion_eta2_hcal->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_cpion_eta3_resolution_fhcal_1oE,Form("2.5<#eta<3.0: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_cpion_eta3_hcal->GetParameter(1),fit_cpion_eta3_hcal->GetParameter(0)),"pl");
  legendNonLin->AddEntry(h_cpion_eta4_resolution_fhcal_1oE,Form("3.0<#eta<4.0: #frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_cpion_eta4_hcal->GetParameter(1),fit_cpion_eta4_hcal->GetParameter(0)),"pl");
  // legendNonLin->AddEntry((TObject*)0,"","");
  legendNonLin->Draw();
  drawLatexAdd("#pi^{-} in FHCAL",0.90,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  cReso2->Print(Form("%s/HCAL_resolution_eta_1oE_Fitted.%s", outputDir.Data(), suffix.Data()));

}


void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value, const Float_t min_towers, const bool normalize, Float_t etalow, Float_t etahigh, Float_t calibconst = 1.0)
{
	Float_t measured_energy;
//      Float_t true_energy;
	t->SetBranchAddress("e", &measured_energy);
//      t->SetBranchAddress("ge", &true_energy);
	Float_t true_eta;
	t->SetBranchAddress("geta", &true_eta);
	Float_t clusterID;
	t->SetBranchAddress("clusterID", &clusterID);
	Float_t nTowers;
	t->SetBranchAddress("ntowers", &nTowers);
	Float_t eventID;
	t->SetBranchAddress("event", &eventID);

	Int_t nentries = Int_t(t->GetEntries());
  Float_t lastEvt = 0;
	for (Int_t i = 0; i < nentries; ++i) {
		if (t->LoadTree(i) < 0)
			break;

		t->GetEntry(i);
    // if(eventID != lastEvt)cout << eventID << "\tE/Etrue: " << measured_energy/true_energy << endl;
    // else cout << "\tE/Etrue: " << measured_energy/true_energy << endl;
		if (
      // ((true_eta > -0.5 && true_eta < 0.5)
      // || (true_eta > -3 && true_eta < -2)
      // || 
      (true_eta > etalow && true_eta < etahigh) //)  //1.7 < 3
    && (measured_energy > min_value && true_energy > 0.1) // above MIP
    // && clusterID==1
    && nTowers>=min_towers // actual cluster with >1 tower
    ){
			h->Fill(measured_energy / calibconst / true_energy);
    }
    // if(eventID != lastEvt)lastEvt = eventID;

	}
	if (normalize)
		h->Scale(1 / h->GetEntries());

	h->SetXTitle("E_{cluster} / E_{true}");
	h->SetYTitle("entries / #scale[0.5]{#sum} entries      ");
}

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs)
{
  if(numInputs<2){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 1, gStyle->GetCanvasDefH()*1);
	  // cPNG->Divide(2, 2);
  } else if(numInputs<5){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 2, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(2, 2);
  }else if(numInputs<7){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(3, 2);
  }else if(numInputs<10){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(3, 3);
  } else if(numInputs<13){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(4, 3);
  } else if(numInputs<17){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(4, 4);
  } else if(numInputs<21){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(5, 4);
  } else {
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*5);
	  cPNG->Divide(5, 5);
  }

}
