#include "plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value,const Float_t min_towers, const bool normalize, Float_t etalow, Float_t etahigh);
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);

void studies_leak2(
  Double_t jetradiusplot = 0.7,
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "e+p: 10#times250 GeV^{2}";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TString inputFolderNames[4] = {"Fun4All_G4_FullDetector_Pythia_nobeampipe","Fun4All_G4_FullDetector_Pythia_beampipe","Fun4All_G4_HighAccDetector_Pythia","Fun4All_G4_FullEtaLeakTestDetector_Pythia"};
  TString namePlotInput[4] = {"w/o pipe","w/ pipe","#eta_{max}=4","#eta_{max}=5.1"};


  //************************** Read data **************************************************
  const Int_t nEne = 23;
  const static Double_t partE[]   = {3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0,  10., 11., 13., 15., 20., 25., 30., 40., 50., 60., 70., 80., 150.};
  Int_t useE[nEne][4]                = {{ 0}};
  Double_t fitrange_low[nEne]     = {0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6,0.65, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7};
  Double_t fitrange_hig[nEne]     = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0};
  TH1F* h_pion_fhcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta1[nEne][4] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta2[nEne][4] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta3[nEne][4] 	    = {NULL};
  TH1F* h_pion_fhcal_E_eta4[nEne][4] 	    = {NULL};
  TH1F* h_pion_fhcal_E_allclus[nEne][4] 	    = {NULL};
  TH1F* h_gamma_fhcal_E[nEne][4] 	    = {NULL};
  Int_t activeE[4] = {0};
  Color_t colorInput[4] = {kGreen+2,kOrange+2,kBlue+2,kRed+2};
  Marker_t markerInput[4] = {27,28,46,20};
  TH1F* h_pion__dEta_hcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dPhi_hcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dR_hcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dEta_ecal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dPhi_ecal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dR_ecal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dEta_ecalfhcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dPhi_ecalfhcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion__dR_ecalfhcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion_ecalfhcal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion_ecal_E[nEne][4] 	    = {NULL};
  TH1F* h_pion_ecalfhcalcomb_E[nEne][4] 	    = {NULL};
  TH2F* h2D_pion_ecalfhcal_E[nEne][4] 	    = {NULL};
  TH2F* h2D_NTowers_HCAL_pi[nEne][4] 	    = {NULL};
  TH2F* h2D_NTowers_ECAL_pi[nEne][4] 	    = {NULL};

  Double_t etalowdef = 3.0;
  Double_t etahighdef = 4.0;

  for(Int_t iInp=0; iInp<4;iInp++){
    for(Int_t epart=0; epart<NELEMS(partE); epart++){
      if(partE[epart]!=15) continue;
      cout << partE[epart] << endl;
      // TTree *const t_pion_fhcal =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),partE[epart]), "READ"))->Get("ntp_cluster"); //ntp_cluster
      TChain* t_pion_fhcal = new TChain("ntp_gshower", "ntp_gshower");
      t_pion_fhcal->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),partE[epart]));
      for(Int_t ii=1; ii<4;ii++){
          t_pion_fhcal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/800_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      }
      for(Int_t ii=2; ii<4;ii++){
          t_pion_fhcal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      }
      if(t_pion_fhcal){
        useE[epart][iInp] = 1;
        activeE[iInp]++;
        h_pion_fhcal_E[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta1[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta1_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta2[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta2_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta3[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta3_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta4[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta4_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_allclus[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_allclus_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_fhcal_E[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
        fill_histogram(h_pion_fhcal_E_eta1[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 1, 1.7);
        fill_histogram(h_pion_fhcal_E_eta2[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 1.7, 2.5);
        fill_histogram(h_pion_fhcal_E_eta3[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 2.5, 3.0);
        fill_histogram(h_pion_fhcal_E_eta4[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 3.0, 4.0);
        fill_histogram(h_pion_fhcal_E_allclus[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 1, true, etalowdef, etahighdef);
      }
      TTree *const t_pion_fecalfhcal =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),partE[epart]), "READ"))->Get("ntp_cluster");
      TTree *const t_pion_fecal =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),partE[epart]), "READ"))->Get("ntp_cluster");
      TTree *const t_pion_fecalfhcal_hcal =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),partE[epart]), "READ"))->Get("ntp_cluster");
      TTree *const t_pion_fecalfhcal_ecal =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),partE[epart]), "READ"))->Get("ntp_cluster");
      if(t_pion_fecalfhcal){
        h_pion_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_ecalfhcal_E[epart][iInp], t_pion_fecalfhcal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_pion_fecal){
        h_pion_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion_ecal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_ecal_E[epart][iInp], t_pion_fecal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_pion_fecalfhcal_hcal && t_pion_fecalfhcal_ecal){
        h_pion__dEta_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dEta_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dPhi_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dR_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dEta_hcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dEta_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_hcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dPhi_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_hcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dR_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dEta_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion__dEta_ecal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion__dPhi_ecal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion__dR_ecal_E_%d",epart), "", 200, -0.5, 0.5);


        h2D_NTowers_HCAL_pi[epart][iInp] 	= new TH2F(Form("h2D_NTowers_HCAL_pi%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_NTowers_ECAL_pi[epart][iInp] 	= new TH2F(Form("h2D_NTowers_ECAL_pi%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_pion_ecalfhcal_E[epart][iInp] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0,partE[epart]*1.15, 100, 0,partE[epart]*1.15);
        h_pion_ecalfhcalcomb_E[epart][iInp] 	= new TH1F(Form("h_pion_ecalfhcalcomb_E_%d",epart), "", 100, 0.0, 2.0);
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
          if (((true_eta1 > 1.7 && true_eta1 < 3.0)
            || (true_eta2 > 1.7 && true_eta2 < 3.0))  //1.7 < 3
          && (measured_energy1 > 0.1) // above MIP
          && (measured_energy1 > 0.1) // above MIP
          // && nTowers1>=2 && nTowers2>=2 // actual cluster with >1 tower
          ){
            h2D_NTowers_HCAL_pi[epart][iInp]->Fill(measured_energy1, nTowers1);
            h2D_NTowers_ECAL_pi[epart][iInp]->Fill(measured_energy2, nTowers2);
            h2D_pion_ecalfhcal_E[epart][iInp]->Fill(measured_energy2, measured_energy1);
            h_pion_ecalfhcalcomb_E[epart][iInp]->Fill((measured_energy2+measured_energy1)/partE[epart]);
            h_pion__dEta_ecalfhcal_E[epart][iInp]->Fill(rec_eta1-rec_eta2);
            h_pion__dPhi_ecalfhcal_E[epart][iInp]->Fill(rec_phi1-rec_phi2);
            h_pion__dR_ecalfhcal_E[epart][iInp]->Fill(TMath::Sqrt(TMath::Power(rec_phi1-rec_phi2,2)+TMath::Power(rec_eta1-rec_eta2,2)));
            h_pion__dEta_hcal_E[epart][iInp]->Fill(rec_eta1-true_eta1);
            h_pion__dPhi_hcal_E[epart][iInp]->Fill(rec_phi1-true_phi1);
            h_pion__dR_hcal_E[epart][iInp]->Fill(TMath::Sqrt(TMath::Power(rec_phi1-true_phi1,2)+TMath::Power(rec_eta1-true_eta1,2)));
            h_pion__dEta_ecal_E[epart][iInp]->Fill(rec_eta2-true_eta2);
            h_pion__dPhi_ecal_E[epart][iInp]->Fill(rec_phi2-true_phi2);
            h_pion__dR_ecal_E[epart][iInp]->Fill(TMath::Sqrt(TMath::Power(rec_phi2-true_phi2,2)+TMath::Power(rec_eta2-true_eta2,2)));
          }
        }
        h_pion_ecalfhcalcomb_E[epart][iInp]->Scale(1 / h_pion_ecalfhcalcomb_E[epart][iInp]->GetEntries());
        h_pion__dEta_ecalfhcal_E[epart][iInp]->Scale(1 / h_pion__dEta_ecalfhcal_E[epart][iInp]->GetEntries());
        h_pion__dPhi_ecalfhcal_E[epart][iInp]->Scale(1 / h_pion__dPhi_ecalfhcal_E[epart][iInp]->GetEntries());
        h_pion__dR_ecalfhcal_E[epart][iInp]->Scale(1 / h_pion__dR_ecalfhcal_E[epart][iInp]->GetEntries());
        h_pion__dEta_hcal_E[epart][iInp]->Scale(1 / h_pion__dEta_hcal_E[epart][iInp]->GetEntries());
        h_pion__dPhi_hcal_E[epart][iInp]->Scale(1 / h_pion__dPhi_hcal_E[epart][iInp]->GetEntries());
        h_pion__dR_hcal_E[epart][iInp]->Scale(1 / h_pion__dR_hcal_E[epart][iInp]->GetEntries());
        h_pion__dEta_ecal_E[epart][iInp]->Scale(1 / h_pion__dEta_ecal_E[epart][iInp]->GetEntries());
        h_pion__dPhi_ecal_E[epart][iInp]->Scale(1 / h_pion__dPhi_ecal_E[epart][iInp]->GetEntries());
        h_pion__dR_ecal_E[epart][iInp]->Scale(1 / h_pion__dR_ecal_E[epart][iInp]->GetEntries());
      }
    }
  }



  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
	split_canvas(cPNG, "cPNG1", activeE[0]);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 1.5,1000,0, 0.070);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_reso[nEne][4] = {NULL};
  TF1 *fit_reso_NCell[nEne][4] = {NULL};
  TF1 *fit_reso_eta1[nEne][4] = {NULL};
  TF1 *fit_reso_eta2[nEne][4] = {NULL};
  TF1 *fit_reso_eta3[nEne][4] = {NULL};
  TF1 *fit_reso_eta4[nEne][4] = {NULL};

  TH1F*  h_pion_resolution_fhcal_E 	= new TH1F("h_pion_resolution_fhcal_E", "", 400, 0.0, 200.0);
  TH1F*  h_pion_resolution_fhcal_1oE 	= new TH1F("h_pion_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  Int_t padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart][0]) continue;
  	cPNG->SetRightMargin(0.01);
  	cPNG->SetTopMargin(0.01);
  	cPNG->SetLeftMargin(0.08);
  	cPNG->SetBottomMargin(0.1);
  	cPNG->cd(padnum+1);
    histEPlotDummy->DrawCopy();
    drawLatexAdd(Form("%1.1f GeV #pi^{-} in %1.1f<#eta<%1.1f",partE[epart],etalowdef,etahighdef),0.15,0.15,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("rec. w/ FHCAL"),0.95,0.15,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    TLegend* legendHE_gamma           = GetAndSetLegend2(0.30, 0.95-(5*textSizeLabelsRel), 0.95, 0.95,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    for(Int_t iInp=0; iInp<4;iInp++){
      if(!useE[epart][iInp]) continue;

      DrawGammaSetMarker(h_pion_fhcal_E[epart][iInp], markerInput[iInp], 2, colorInput[iInp], colorInput[iInp]);
      h_pion_fhcal_E[epart][iInp]->Draw("same");

      fit_reso[epart][iInp] = new TF1(Form("fit_reso_%d",epart), "gaus", iInp<2 ? fitrange_low[epart]-0.1 : fitrange_low[epart],iInp<2 ? fitrange_hig[epart]-0.1 : fitrange_hig[epart]);

      h_pion_fhcal_E[epart][iInp]->Fit(fit_reso[epart][iInp],"L0RMEQ");
      fit_reso[epart][iInp]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso[epart][iInp], 1, 2, colorInput[iInp]+2);
      fit_reso[epart][iInp]->Draw("same");
      legendHE_gamma->AddEntry(h_pion_fhcal_E[epart][iInp],Form("%s: mean: %1.2f, #sigma/#it{E}: %1.2f%%",namePlotInput[iInp].Data(),fit_reso[epart][iInp]->GetParameter(1),100*fit_reso[epart][iInp]->GetParameter(2)/fit_reso[epart][iInp]->GetParameter(1)),"pl");
    }
    legendHE_gamma->Draw();

    // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // h_pion_resolution_fhcal_E->SetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(2)/fit_reso[epart]->GetParameter(2),2)+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2))); //+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)
    // h_pion_resolution_fhcal_1oE->SetBinContent(h_pion_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_resolution.%s", outputDir.Data(), suffix.Data()));
 


}


void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value, const Float_t min_towers, const bool normalize, Float_t etalow, Float_t etahigh)
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
    // && nTowers>=min_towers // actual cluster with >1 tower
    ){
			h->Fill(measured_energy / true_energy);
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
  }else if(numInputs<5){
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


/*

  Double_t etalowdef = 3.0;
  Double_t etahighdef = 4.0;

  for(Int_t iInp=0; iInp<4;iInp++){
    for(Int_t epart=0; epart<NELEMS(partE); epart++){
      cout << partE[epart] << endl;
      // TChain* t_pion_fhcal = new TChain("ntp_cluster", "ntp_cluster");
      TChain* t_pion_fhcal = new TChain("ntp_gshower", "ntp_gshower");
      t_pion_fhcal->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),partE[epart]));
      for(Int_t ii=1; ii<4;ii++){
          t_pion_fhcal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/800_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      }
      for(Int_t ii=2; ii<4;ii++){
          t_pion_fhcal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      }
      if(t_pion_fhcal){
        useE[epart][iInp] = 1;
        activeE[iInp]++;
        h_pion_fhcal_E[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta1[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta1_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta2[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta2_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta3[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta3_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_eta4[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_eta4_%d",epart), "", 100, 0.0, 2.0);
        h_pion_fhcal_E_allclus[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_allclus_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_fhcal_E[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
        fill_histogram(h_pion_fhcal_E_eta1[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 1, 1.7);
        fill_histogram(h_pion_fhcal_E_eta2[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 1.7, 2.5);
        fill_histogram(h_pion_fhcal_E_eta3[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 2.5, 3.0);
        fill_histogram(h_pion_fhcal_E_eta4[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 2, true, 3.0, 4.0);
        fill_histogram(h_pion_fhcal_E_allclus[epart][iInp], t_pion_fhcal, partE[epart], 0.3, 1, true, etalowdef, etahighdef);
      }
      TChain* t_pion_fecalfhcal = new TChain("ntp_cluster", "ntp_cluster");
      TChain* t_pion_fecal = new TChain("ntp_cluster", "ntp_cluster");
      TChain* t_pion_fecalfhcal_hcal = new TChain("ntp_cluster", "ntp_cluster");
      TChain* t_pion_fecalfhcal_ecal = new TChain("ntp_cluster", "ntp_cluster");
      t_pion_fecalfhcal->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),partE[epart]));
      t_pion_fecal->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),partE[epart]));
      t_pion_fecalfhcal_hcal->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),partE[epart]));
      t_pion_fecalfhcal_ecal->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),partE[epart]));

      // for(Int_t ii=2; ii<4;ii++){
      //   // if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/800_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]))){
      //     t_pion_fecalfhcal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      //     t_pion_fecalfhcal_hcal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      //   // }
      //   // if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/800_%d_%1.1f/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]))){
      //     t_pion_fecal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      //     t_pion_fecalfhcal_ecal->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4femc_eval.root",inputFolderNames[iInp].Data(),ii,partE[epart]));
      //   // }
      // }

      if(t_pion_fecalfhcal){
        h_pion_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion_ecalfhcal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_ecalfhcal_E[epart][iInp], t_pion_fecalfhcal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_pion_fecal){
        h_pion_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion_ecal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_ecal_E[epart][iInp], t_pion_fecal, partE[epart], 0.3, 2, true, etalowdef, etahighdef);
      }
      if(t_pion_fecalfhcal_hcal && t_pion_fecalfhcal_ecal){
        h_pion__dEta_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dEta_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dPhi_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_ecalfhcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dR_ecalfhcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dEta_hcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dEta_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_hcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dPhi_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_hcal_E[epart][iInp] 	= new TH1F(Form("h_pion__dR_hcal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dEta_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion__dEta_ecal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dPhi_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion__dPhi_ecal_E_%d",epart), "", 200, -0.5, 0.5);
        h_pion__dR_ecal_E[epart][iInp] 	= new TH1F(Form("h_pion__dR_ecal_E_%d",epart), "", 200, -0.5, 0.5);


        h2D_NTowers_HCAL_pi[epart][iInp] 	= new TH2F(Form("h2D_NTowers_HCAL_pi%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_NTowers_ECAL_pi[epart][iInp] 	= new TH2F(Form("h2D_NTowers_ECAL_pi%d",epart), "", 100, 0,partE[epart]*1.15, 40, -0.5,39.5);
        h2D_pion_ecalfhcal_E[epart][iInp] 	= new TH2F(Form("h2D_pion_ecalfhcal_E_%d",epart), "", 100, 0,partE[epart]*1.15, 100, 0,partE[epart]*1.15);
        h_pion_ecalfhcalcomb_E[epart][iInp] 	= new TH1F(Form("h_pion_ecalfhcalcomb_E_%d",epart), "", 100, 0.0, 2.0);
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
          if (((true_eta1 > 1.7 && true_eta1 < 3.0)
            || (true_eta2 > 1.7 && true_eta2 < 3.0))  //1.7 < 3
          && (measured_energy1 > 0.1) // above MIP
          && (measured_energy1 > 0.1) // above MIP
          // && nTowers1>=2 && nTowers2>=2 // actual cluster with >1 tower
          ){
            h2D_NTowers_HCAL_pi[epart][iInp]->Fill(measured_energy1, nTowers1);
            h2D_NTowers_ECAL_pi[epart][iInp]->Fill(measured_energy2, nTowers2);
            h2D_pion_ecalfhcal_E[epart][iInp]->Fill(measured_energy2, measured_energy1);
            h_pion_ecalfhcalcomb_E[epart][iInp]->Fill((measured_energy2+measured_energy1)/partE[epart]);
            h_pion__dEta_ecalfhcal_E[epart][iInp]->Fill(rec_eta1-rec_eta2);
            h_pion__dPhi_ecalfhcal_E[epart][iInp]->Fill(rec_phi1-rec_phi2);
            h_pion__dR_ecalfhcal_E[epart][iInp]->Fill(TMath::Sqrt(TMath::Power(rec_phi1-rec_phi2,2)+TMath::Power(rec_eta1-rec_eta2,2)));
            h_pion__dEta_hcal_E[epart][iInp]->Fill(rec_eta1-true_eta1);
            h_pion__dPhi_hcal_E[epart][iInp]->Fill(rec_phi1-true_phi1);
            h_pion__dR_hcal_E[epart][iInp]->Fill(TMath::Sqrt(TMath::Power(rec_phi1-true_phi1,2)+TMath::Power(rec_eta1-true_eta1,2)));
            h_pion__dEta_ecal_E[epart][iInp]->Fill(rec_eta2-true_eta2);
            h_pion__dPhi_ecal_E[epart][iInp]->Fill(rec_phi2-true_phi2);
            h_pion__dR_ecal_E[epart][iInp]->Fill(TMath::Sqrt(TMath::Power(rec_phi2-true_phi2,2)+TMath::Power(rec_eta2-true_eta2,2)));
          }
        }
        h_pion_ecalfhcalcomb_E[epart][iInp]->Scale(1 / h_pion_ecalfhcalcomb_E[epart][iInp]->GetEntries());
        h_pion__dEta_ecalfhcal_E[epart][iInp]->Scale(1 / h_pion__dEta_ecalfhcal_E[epart][iInp]->GetEntries());
        h_pion__dPhi_ecalfhcal_E[epart][iInp]->Scale(1 / h_pion__dPhi_ecalfhcal_E[epart][iInp]->GetEntries());
        h_pion__dR_ecalfhcal_E[epart][iInp]->Scale(1 / h_pion__dR_ecalfhcal_E[epart][iInp]->GetEntries());
        h_pion__dEta_hcal_E[epart][iInp]->Scale(1 / h_pion__dEta_hcal_E[epart][iInp]->GetEntries());
        h_pion__dPhi_hcal_E[epart][iInp]->Scale(1 / h_pion__dPhi_hcal_E[epart][iInp]->GetEntries());
        h_pion__dR_hcal_E[epart][iInp]->Scale(1 / h_pion__dR_hcal_E[epart][iInp]->GetEntries());
        h_pion__dEta_ecal_E[epart][iInp]->Scale(1 / h_pion__dEta_ecal_E[epart][iInp]->GetEntries());
        h_pion__dPhi_ecal_E[epart][iInp]->Scale(1 / h_pion__dPhi_ecal_E[epart][iInp]->GetEntries());
        h_pion__dR_ecal_E[epart][iInp]->Scale(1 / h_pion__dR_ecal_E[epart][iInp]->GetEntries());
      }
    }
  }



  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
	split_canvas(cPNG, "cPNG1", activeE[0]);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 1.5,1000,0, 0.070);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_reso[nEne][4] = {NULL};
  TF1 *fit_reso_NCell[nEne][4] = {NULL};
  TF1 *fit_reso_eta1[nEne][4] = {NULL};
  TF1 *fit_reso_eta2[nEne][4] = {NULL};
  TF1 *fit_reso_eta3[nEne][4] = {NULL};
  TF1 *fit_reso_eta4[nEne][4] = {NULL};

  TH1F*  h_pion_resolution_fhcal_E 	= new TH1F("h_pion_resolution_fhcal_E", "", 400, 0.0, 200.0);
  TH1F*  h_pion_resolution_fhcal_1oE 	= new TH1F("h_pion_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  Int_t padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart][0]) continue;
  	cPNG->cd(padnum+1);
    histEPlotDummy->DrawCopy();
    drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    for(Int_t iInp=0; iInp<4;iInp++){
      if(!useE[epart][iInp]) continue;

      DrawGammaSetMarker(h_pion_fhcal_E[epart][iInp], markerInput[iInp], 2, colorInput[iInp], colorInput[iInp]);
      h_pion_fhcal_E[epart][iInp]->Draw("same");

      fit_reso[epart][iInp] = new TF1(Form("fit_reso_%d",epart), "gaus", fitrange_low[epart],fitrange_hig[epart]);

      h_pion_fhcal_E[epart][iInp]->Fit(fit_reso[epart][iInp],"L0RMEQ");
      fit_reso[epart][iInp]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso[epart][iInp], 1, 2, colorInput[iInp]+2);
      fit_reso[epart][iInp]->Draw("same");
    }

    // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // h_pion_resolution_fhcal_E->SetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(2)/fit_reso[epart]->GetParameter(2),2)+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2))); //+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)
    // h_pion_resolution_fhcal_1oE->SetBinContent(h_pion_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    padnum++;
  }
  cPNG->Print(Form("%s/HCAL_resolution.%s", outputDir.Data(), suffix.Data()));

*/