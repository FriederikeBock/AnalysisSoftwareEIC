#include "../common/plottingheader.h"
#include "../common/binningheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value,const Float_t min_towers, const bool normalize);


void calibration(
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "e-p, 20 x 250 GeV^{2}";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("calibrationPlots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDir+"/EtaBins");

  TString str_calorimeter[5] = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC"};
  const int _active_calo = 2;
  enum clusterizertype {
      kV1         = 0,
      kV3         = 1,
      k3x3        = 2,
      k5x5        = 3,
      kC3         = 4,
      kC5         = 5,
      kMA         = 6,
      kDummy      = 7
  };
  TString str_clusterizer[7] = {"V1", "V3", "3x3", "5x5", "C3", "C5", "MA"};
  const int _active_algo = 7;
  int marker_clusterizer[_active_algo]     = {20, 22, 29, 47, 33, 20, 22};
  int linestyle_clusterizer[_active_algo]     = {2, 4, 8, 9, 5, 4, 8};


  TFile* inputFileCALOGamma = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CALOSTANDALONE_SimplePhoton/output_RH.root");
  TFile* inputFileCALOPion = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CALOSTANDALONE_SimplePion/output_RH.root");
  TFile* inputFileFHCALPion = new TFile("/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/FHCALSTANDALONE_SimplePion/output_RH.root");
  cout << "using inputs from treeProcessing_calibration!!!!!" << endl;

  TH2F*  h_RH_Reso_gamma_E[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_gamma_E_Eta[_active_calo][_active_algo][nEta+1];
  TH2F*  h_RH_Reso_gamma_Erec[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_gammaCalib_E[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_gammaCalib_Erec[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pion_E[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pion_smear_E[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pion_E_Eta[_active_calo][_active_algo][nEta+1];
  TH2F*  h_RH_Reso_pion_Erec[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pionComb_E[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pionComb_Erec[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pionCombCalib_E[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pionCombCalib_Erec[_active_calo][_active_algo];
  TH2F*  h_RH_Reso_pionCombCalib_smeared_Erec[_active_calo][_active_algo];


  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      for (Int_t eT = 0; eT < nEta+1; eT++){
        h_RH_Reso_gamma_E_Eta[icalo][ialgo][eT] 	= (TH2F*) inputFileCALOGamma->Get(Form("h_RH_Reso_gamma_E_Eta_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),eT));
        h_RH_Reso_pion_E_Eta[icalo][ialgo][eT] 	= (TH2F*) inputFileFHCALPion->Get(Form("h_RH_Reso_pion_E_Eta_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),eT));
      }
      h_RH_Reso_gamma_E[icalo][ialgo] 	= (TH2F*) inputFileCALOGamma->Get(Form("h_RH_Reso_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_gamma_Erec[icalo][ialgo] 	= (TH2F*) inputFileCALOGamma->Get(Form("h_RH_Reso_gamma_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_gammaCalib_E[icalo][ialgo] 	= (TH2F*) inputFileCALOGamma->Get(Form("h_RH_ResoCalib_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_gammaCalib_Erec[icalo][ialgo] 	= (TH2F*) inputFileCALOGamma->Get(Form("h_RH_ResoCalib_gamma_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));

      h_RH_Reso_pion_E[icalo][ialgo] 	= (TH2F*) inputFileFHCALPion->Get(Form("h_RH_Reso_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      // h_RH_Reso_pion_smear_E[icalo][ialgo] 	= (TH2F*) inputFileFHCALPion->Get(Form("h_RH_Reso_pion_smear_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));

      h_RH_Reso_pion_Erec[icalo][ialgo] 	= (TH2F*) inputFileFHCALPion->Get(Form("h_RH_Reso_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_pionComb_E[icalo][ialgo] 	= (TH2F*) inputFileCALOPion->Get(Form("h_RH_Reso_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_pionComb_Erec[icalo][ialgo] 	= (TH2F*) inputFileCALOPion->Get(Form("h_RH_Reso_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_pionCombCalib_E[icalo][ialgo] 	= (TH2F*) inputFileCALOPion->Get(Form("h_RH_ResoCalib_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_pionCombCalib_Erec[icalo][ialgo] 	= (TH2F*) inputFileCALOPion->Get(Form("h_RH_ResoCalib_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
      h_RH_Reso_pionCombCalib_smeared_Erec[icalo][ialgo] 	= (TH2F*) inputFileCALOPion->Get(Form("h_RH_ResoCalib_smeared_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()));
    }
  }

  //************************** Read data **************************************************
  const Int_t nEne = 23;
  const static float partE[]   = {5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0, 10., 11., 13.,      15., 20., 25., 30., 40., 50., 60., 70., 80., 100,      125, 150, 200};
  TH1F*  h_RH_Reso_gamma_bin_E[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_gamma_bin_Erec[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_gammaCalib_bin_E[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_gammaCalib_bin_Erec[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_pion_bin_E[_active_calo][_active_algo][nEne]= {NULL};
  // TH1F*  h_RH_Reso_pion_bin_smear_E[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_pion_bin_E_Eta[_active_calo][_active_algo][nEta+1][nEne]= {NULL};
  TH1F*  h_RH_Reso_pion_bin_Erec[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_pionComb_bin_E[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_pionComb_bin_Erec[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_pionCombCalib_bin_E[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_pionCombCalib_bin_Erec[_active_calo][_active_algo][nEne]= {NULL};
  TH1F*  h_RH_Reso_pionCombCalib_bin_smeared_Erec[_active_calo][_active_algo][nEne]= {NULL};

  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
          Int_t startBin  = h_RH_Reso_gamma_E[icalo][ialgo] ->GetXaxis()->FindBin(partE[epart]+0.001);
          Int_t endBin    = h_RH_Reso_gamma_E[icalo][ialgo] ->GetXaxis()->FindBin(partE[epart+1]-0.001);

          h_RH_Reso_gamma_bin_E[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_gamma_E[icalo][ialgo] ->ProjectionY(Form("h_gamma_E_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_gamma_bin_E[icalo][ialgo][epart]->Scale(1./h_RH_Reso_gamma_bin_E[icalo][ialgo][epart]->GetEntries());
          h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_gamma_Erec[icalo][ialgo] ->ProjectionY(Form("h_gamma_Erec_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]->Scale(1./h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]->GetEntries());

          h_RH_Reso_gammaCalib_bin_E[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_gammaCalib_E[icalo][ialgo] ->ProjectionY(Form("h_gammaCalib_E_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_gammaCalib_bin_E[icalo][ialgo][epart]->Scale(1./h_RH_Reso_gammaCalib_bin_E[icalo][ialgo][epart]->GetEntries());
          h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_gammaCalib_Erec[icalo][ialgo] ->ProjectionY(Form("h_gammaCalib_Erec_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]->Scale(1./h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]->GetEntries());

          for (Int_t eT = 0; eT < nEta+1; eT++){
            h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart] = (TH1F*)h_RH_Reso_pion_E_Eta[icalo][ialgo][eT] ->ProjectionY(Form("h_pion_E_Eta_%s_%s_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),eT, epart),startBin,endBin,"e");
            h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart]->Scale(1./h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart]->GetEntries());
          }
          h_RH_Reso_pion_bin_E[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pion_E[icalo][ialgo] ->ProjectionY(Form("h_pion_E_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_pion_bin_E[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pion_bin_E[icalo][ialgo][epart]->GetEntries());

          // h_RH_Reso_pion_bin_smear_E[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pion_smear_E[icalo][ialgo] ->ProjectionY(Form("h_pion_E_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          // h_RH_Reso_pion_bin_smear_E[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pion_bin_smear_E[icalo][ialgo][epart]->GetEntries());

          h_RH_Reso_pion_bin_Erec[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pion_Erec[icalo][ialgo] ->ProjectionY(Form("h_pion_Erec_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_pion_bin_Erec[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pion_bin_Erec[icalo][ialgo][epart]->GetEntries());

          h_RH_Reso_pionComb_bin_E[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pionComb_E[icalo][ialgo] ->ProjectionY(Form("h_pionComb_E_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_pionComb_bin_E[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pionComb_bin_E[icalo][ialgo][epart]->GetEntries());
          h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pionComb_Erec[icalo][ialgo] ->ProjectionY(Form("h_pionComb_Erec_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]->GetEntries());

          h_RH_Reso_pionCombCalib_bin_E[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pionCombCalib_E[icalo][ialgo] ->ProjectionY(Form("h_pionCombCalib_E_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_pionCombCalib_bin_E[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pionCombCalib_bin_E[icalo][ialgo][epart]->GetEntries());
          h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pionCombCalib_Erec[icalo][ialgo] ->ProjectionY(Form("h_RH_Reso_pionCombCalib_bin_Erec_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
          h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]->GetEntries());
          if(h_RH_Reso_pionCombCalib_smeared_Erec[icalo][ialgo]){
            h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart] = (TH1F*)h_RH_Reso_pionCombCalib_smeared_Erec[icalo][ialgo] ->ProjectionY(Form("h_RH_Reso_pionCombCalib_bin_smeared_Erec_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),startBin,endBin,"e");
            h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart]->Scale(1./h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart]->GetEntries());
          }

      }
    }
  }

  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
	split_canvas(cPNG, "cPNG1", nEne);
  // cPNG->SetLogy();
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  // TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0.2, 1.5,1000,0, 0.07);
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0.2, 1.5,1000,0.0001, 0.06);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);


  TF1 *fit_reso_pion_Eta[_active_algo][nEta+1][nEne] = {NULL};
  TH1F*  h_FHCALSTANDALONE_pion_reso_1oE_Eta[_active_calo][_active_algo][nEta+1] 	= {NULL};
  TH1F*  h_FHCALSTANDALONE_pion_meanE_E_Eta[_active_calo][_active_algo][nEta+1] 	= {NULL};
  for (Int_t eT = 9; eT < 14; eT++){
    for(int icalo=0;icalo<_active_calo;icalo++){
      if(icalo!=0)continue;
      for(int ialgo=0;ialgo<_active_algo;ialgo++){
        if(ialgo==kC3 || ialgo==k3x3) continue;
        if(!(ialgo==kMA)) continue;
        Int_t padnum =0;
        for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
          DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
          cPNG->cd(padnum+1);

          histEPlotDummy->DrawCopy();
          if(h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart]){
            float scalerange = 1;
            DrawGammaSetMarker(h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart], marker_clusterizer[ialgo], 3, colorClus[1][ialgo], colorClus[1][ialgo]);
            h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart]->Draw("same");
            fit_reso_pion_Eta[ialgo][eT][epart]    = new TF1(Form("fit_reso_pion_Eta_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                                      0.3, 1.0);
            fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(0, 0.05);
            fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(1, 0.73);
            fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(2, 0.01);
            fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(3, 0.015);
            h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart]->Fit(fit_reso_pion_Eta[ialgo][eT][epart],"L0MEQ", "", 0.35, 0.95);
            if(fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(1)>1 || fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(1)< 0.5){
              fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(0, 0.05);
              fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(1, 0.73);
              fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(2, 0.01);
              fit_reso_pion_Eta[ialgo][eT][epart]->SetParameter(3, 0.015);
              h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart]->Fit(fit_reso_pion_Eta[ialgo][eT][epart],"L0MEQ", "", 0.3, 1.0);
            }
            if(fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(1)<0 || fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(2)<0.015){
              if(ialgo==kV1)
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", 0.15+0.04*TMath::Log(partE[epart]),0.45);
              else if(ialgo==kV3)
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", scalerange*(0.15+0.1*TMath::Log(partE[epart])),scalerange*(0.85));
              else if(ialgo==k3x3)
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.80));
              else if(ialgo==k5x5)
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
              else if(ialgo==kC3)
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.75));
              else if(ialgo==kC5)
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
              else if(ialgo==kMA)
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", scalerange*(0.35+0.145*TMath::Log(TMath::Sqrt(partE[epart]))),scalerange*(0.85));
              else
                fit_reso_pion_Eta[ialgo][eT][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,eT,epart), "gaus", scalerange*(0.25+0.1*TMath::Log(partE[epart])),scalerange*(0.9));

              h_RH_Reso_pion_bin_E_Eta[icalo][ialgo][eT][epart]->Fit(fit_reso_pion_Eta[ialgo][eT][epart],"L0RMEQ");
            }
            fit_reso_pion_Eta[ialgo][eT][epart]->SetNpx(10000);
            DrawGammaSetMarkerTF1(fit_reso_pion_Eta[ialgo][eT][epart], linestyle_clusterizer[ialgo], 2, colorClus[1][ialgo]+2);
            fit_reso_pion_Eta[ialgo][eT][epart]->Draw("same");
            if(!h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT])h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT] 	= new TH1F(Form("h_pion_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
            if(!(eT==13 && partE[epart]<15)){
              if(!(eT==12 && (partE[epart]==7 || partE[epart]==10))){
              h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->SetBinContent(h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(2)/fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(1));
              }
            }
            if(!h_FHCALSTANDALONE_pion_meanE_E_Eta[icalo][ialgo][eT])h_FHCALSTANDALONE_pion_meanE_E_Eta[icalo][ialgo][eT] 	= new TH1F(Form("h_FHCALSTANDALONE_pion_meanE_E_Eta_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
            h_FHCALSTANDALONE_pion_meanE_E_Eta[icalo][ialgo][eT]->SetBinContent(h_FHCALSTANDALONE_pion_meanE_E_Eta[icalo][ialgo][eT]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(1));

            drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

            drawLatexAdd(Form("mean: %1.2f ",fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(1)),0.15,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
            drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(2)),0.15,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
            drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(2)/fit_reso_pion_Eta[ialgo][eT][epart]->GetParameter(1)),0.15,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          }
          padnum++;
        }
        cPNG->Print(Form("%s/EtaBins/FHCALSTANDALONE_pion_bins_%s%s_Etrue_Eta%d.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),eT, suffix.Data()));
      }
    }
  }

TCanvas* cResoEta2 = new TCanvas("cResoEta2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cResoEta2, 0.095, 0.01, 0.01, 0.105);
  TH2F* histResoEtaDummy                           = new TH2F("histResoEtaDummy","histResoEtaDummy",1000,0, 0.49,1000,0, 79);
  SetStyleHistoTH2ForGraphs(histResoEtaDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoEtaDummy->GetXaxis()->SetNoExponent();
  histResoEtaDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoEtaDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
    histResoEtaDummy->Draw();
    TLegend* legendHE_pion           = GetAndSetLegend2(0.15, 0.95-2*textSizeLabelsRel, 0.7, 0.95-(9*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

    TF1* fit_reso_clusdep_pion_1oE_Eta[nEta+1][_active_algo];
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!(ialgo==kMA)) continue;
      for (Int_t eT = 9; eT < 14; eT++){
        DrawGammaSetMarker(h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT], markerStyleEta[eT], 3, colorEta[eT], colorEta[eT]);
        h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->SetLineColor(colorEta[eT]);
        h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->SetLineStyle(1);
        h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->SetLineWidth(3);
        h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->Draw("same,p");
        fit_reso_clusdep_pion_1oE_Eta[eT][ialgo] = new TF1(Form("fit_reso_clusdep_pion_1oE_%d%d",ialgo,icalo), "[0]+[1]*x+[2]*TMath::Power(x,2)", 0.05,0.47);
        fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->SetParLimits(0,0,30);
        fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->SetParameter(2,0.3);
        fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->FixParameter(2,0.08);
        if(eT==9)
          h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->Fit(fit_reso_clusdep_pion_1oE_Eta[eT][ialgo],"QMNE","",0.1,0.45);
        else if(eT==13)
          h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->Fit(fit_reso_clusdep_pion_1oE_Eta[eT][ialgo],"QMNE","",0.05,0.25);
        else
          h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT]->Fit(fit_reso_clusdep_pion_1oE_Eta[eT][ialgo],"QMNE","",0.08,0.4);
        // fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->SetNpx(10000);
        DrawGammaSetMarkerTF1(fit_reso_clusdep_pion_1oE_Eta[eT][ialgo], 1, 3, colorEta[eT]);
        fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->SetLineStyle(1);
        fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->Draw("same");
        legendHE_pion->AddEntry(h_FHCALSTANDALONE_pion_reso_1oE_Eta[icalo][ialgo][eT],Form("%1.1f < #it{#eta} < %1.1f: #sigma/#it{E} =  %1.2f/#sqrt{#it{E}} #oplus %1.2f", partEta[eT], partEta[eT+1],fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->GetParameter(1),fit_reso_clusdep_pion_1oE_Eta[eT][ialgo]->GetParameter(0)),"pl");
      }
    }
    legendHE_pion->Draw();
    // drawLatexAdd(Form("#pi^{-} in %s",str_calorimeter[icalo].Data()),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#pi^{-} in %s",str_calorimeter[icalo].Data()),0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("placeholder data"),0.5,0.5,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE,kOrange+2);
    drawLatexAdd("MA clusterizer: E_{seed} = 500 MeV, E_{agg} = 100 MeV",0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);


    cResoEta2->Print(Form("%s/EtaBins/FHCALSTANDALONE_pion_reso_1oE_%s_EtaSlices.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
  }

	split_canvas(cPNG, "cPNG1", nEne);

  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0.2, 1.5,1000,0.0001, 0.02);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  TF1 *fit_reso_pion[_active_algo][nEne] = {NULL};
  TH1F*  h_FHCALSTANDALONE_pion_reso_1oE[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FHCALSTANDALONE_pion_meanE_E[_active_calo][_active_algo] 	= {NULL};
  TF1 *fit_reso_pionComb[_active_algo][nEne] = {NULL};
  TH1F*  h_FHCALSTANDALONE_pionComb_reso_1oE[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FHCALSTANDALONE_pionComb_meanE_E[_active_calo][_active_algo] 	= {NULL};
  TF1 *fit_reso_pionCombCalib[_active_algo][nEne] = {NULL};
  TH1F*  h_FHCALSTANDALONE_pionCombCalib_reso_1oE[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FHCALSTANDALONE_pionCombCalib_meanE_E[_active_calo][_active_algo] 	= {NULL};
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
       if(ialgo==kC3 || ialgo==k3x3) continue;
     Int_t padnum =0;
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1);

        histEPlotDummy->DrawCopy();
        if(h_RH_Reso_pion_bin_E[icalo][ialgo][epart]){
          float scalerange = 1;
          DrawGammaSetMarker(h_RH_Reso_pion_bin_E[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[1][ialgo], colorClus[1][ialgo]);
          h_RH_Reso_pion_bin_E[icalo][ialgo][epart]->Draw("same");
          if(ialgo==kV1)
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 0.15+0.04*TMath::Log(partE[epart]),0.45);
          else if(ialgo==kV3)
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.15+0.1*TMath::Log(partE[epart])),scalerange*(0.85));
          else if(ialgo==k3x3)
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.80));
          else if(ialgo==k5x5)
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC3)
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC5)
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kMA)
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.35+0.145*TMath::Log(TMath::Sqrt(partE[epart]))),scalerange*(0.85));
          else
            fit_reso_pion[ialgo][epart] = new TF1(Form("fit_reso_pion_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.25+0.1*TMath::Log(partE[epart])),scalerange*(0.9));

          h_RH_Reso_pion_bin_E[icalo][ialgo][epart]->Fit(fit_reso_pion[ialgo][epart],"L0RMEQ");
          fit_reso_pion[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pion[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[1][ialgo]+2);
          fit_reso_pion[ialgo][epart]->Draw("same");
          if(!h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo])h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo] 	= new TH1F(Form("h_pion_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_pion[ialgo][epart]->GetParameter(2)/fit_reso_pion[ialgo][epart]->GetParameter(1));
          if(!h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo])h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pion_meanE_E_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_pion[ialgo][epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_pion[ialgo][epart]->GetParameter(1)),0.15,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pion[ialgo][epart]->GetParameter(2)),0.15,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pion[ialgo][epart]->GetParameter(2)/fit_reso_pion[ialgo][epart]->GetParameter(1)),0.15,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_pionComb_bin_E[icalo][ialgo][epart]){
          float scalerange = 1;
          DrawGammaSetMarker(h_RH_Reso_pionComb_bin_E[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[0][ialgo], colorClus[0][ialgo]);
          h_RH_Reso_pionComb_bin_E[icalo][ialgo][epart]->Draw("same");
          if(ialgo==kV1)
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", 0.15+0.04*TMath::Log(partE[epart]),0.45);
          else if(ialgo==kV3)
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.15+0.1*TMath::Log(partE[epart])),scalerange*(0.85));
          else if(ialgo==k3x3)
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.80));
          else if(ialgo==k5x5)
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC3)
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC5)
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kMA)
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.35+0.145*TMath::Log(TMath::Sqrt(partE[epart]))),scalerange*(0.85));
          else
            fit_reso_pionComb[ialgo][epart] = new TF1(Form("fit_reso_pionComb_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.25+0.1*TMath::Log(partE[epart])),scalerange*(0.9));

          h_RH_Reso_pionComb_bin_E[icalo][ialgo][epart]->Fit(fit_reso_pionComb[ialgo][epart],"L0RMEQ");
          fit_reso_pionComb[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pionComb[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[0][ialgo]+2);
          fit_reso_pionComb[ialgo][epart]->Draw("same");
          if(!h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo])h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo] 	= new TH1F(Form("h_pionComb_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_pionComb[ialgo][epart]->GetParameter(2)/fit_reso_pionComb[ialgo][epart]->GetParameter(1));
          if(!h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo])h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pionComb_meanE_E_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_pionComb[ialgo][epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_pionComb[ialgo][epart]->GetParameter(1)),0.15,0.65,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pionComb[ialgo][epart]->GetParameter(2)),0.15,0.60,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pionComb[ialgo][epart]->GetParameter(2)/fit_reso_pionComb[ialgo][epart]->GetParameter(1)),0.15,0.55,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_pionCombCalib_bin_E[icalo][ialgo][epart]){
          float scalerange = 1;
          DrawGammaSetMarker(h_RH_Reso_pionCombCalib_bin_E[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[2][ialgo], colorClus[2][ialgo]);
          h_RH_Reso_pionCombCalib_bin_E[icalo][ialgo][epart]->Draw("same");
          if(ialgo==kV1)
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", 0.15+0.04*TMath::Log(partE[epart]),0.45);
          else if(ialgo==kV3)
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.15+0.1*TMath::Log(partE[epart])),scalerange*(0.85));
          else if(ialgo==k3x3)
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.80));
          else if(ialgo==k5x5)
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC3)
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC5)
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kMA)
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", 1.2*(0.35+0.145*TMath::Log(TMath::Sqrt(partE[epart]))),1.2*(0.85));
          else
            fit_reso_pionCombCalib[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo, epart), "gaus", scalerange*(0.25+0.1*TMath::Log(partE[epart])),scalerange*(0.9));

          h_RH_Reso_pionCombCalib_bin_E[icalo][ialgo][epart]->Fit(fit_reso_pionCombCalib[ialgo][epart],"L0RMEQ");
          fit_reso_pionCombCalib[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pionCombCalib[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[2][ialgo]+2);
          fit_reso_pionCombCalib[ialgo][epart]->Draw("same");
          if(!h_FHCALSTANDALONE_pionCombCalib_reso_1oE[icalo][ialgo])h_FHCALSTANDALONE_pionCombCalib_reso_1oE[icalo][ialgo] 	= new TH1F(Form("h_pionCombCalib_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FHCALSTANDALONE_pionCombCalib_reso_1oE[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionCombCalib_reso_1oE[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_pionCombCalib[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib[ialgo][epart]->GetParameter(1));
          if(!h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo])h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pionCombCalib_meanE_E_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_pionCombCalib[ialgo][epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_pionCombCalib[ialgo][epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pionCombCalib[ialgo][epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pionCombCalib[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib[ialgo][epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        padnum++;
    }
      cPNG->Print(Form("%s/FHCALSTANDALONE_pion_bins_%s%s_Etrue.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }


  TF1 *fit_reso_pion_Erec[_active_algo][nEne] = {NULL};
  TH1F*  h_FHCALSTANDALONE_pion_reso_1oErec[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FHCALSTANDALONE_pion_meanE_Erec[_active_calo][_active_algo] 	= {NULL};
  TF1 *fit_reso_pionComb_Erec[_active_algo][nEne] = {NULL};
  TH1F*  h_FHCALSTANDALONE_pionComb_reso_1oErec[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FHCALSTANDALONE_pionComb_meanE_Erec[_active_calo][_active_algo] 	= {NULL};
  TF1 *fit_reso_pionCombCalib_Erec[_active_algo][nEne] = {NULL};
  TF1 *fit_reso_pionCombCalib_smeared_Erec[_active_algo][nEne] = {NULL};
  TH1F*  h_FHCALSTANDALONE_pionCombCalib_reso_1oErec[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[_active_calo][_active_algo] 	= {NULL};
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      Int_t padnum =0;
      if(ialgo==kC3 || ialgo==k3x3) continue;
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1);

        histEPlotDummy->DrawCopy();
        if(h_RH_Reso_pion_bin_Erec[icalo][ialgo][epart]){
          float scalerange = 1.0;
          DrawGammaSetMarker(h_RH_Reso_pion_bin_Erec[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[1][ialgo], colorClus[1][ialgo]);
          // h_RH_Reso_pion_bin_Erec[icalo][ialgo][epart]->Draw("same");
          if(ialgo==kV1)
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 0.15+0.04*TMath::Log(partE[epart]),0.45);
          else if(ialgo==kV3)
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.15+0.1*TMath::Log(partE[epart])),scalerange*(0.85));
          else if(ialgo==k3x3)
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.80));
          else if(ialgo==k5x5)
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC3)
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC5)
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kMA)
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.35+0.145*TMath::Log(TMath::Sqrt(partE[epart]))),scalerange*(0.85));
          else
            fit_reso_pion_Erec[ialgo][epart] = new TF1(Form("fit_reso_pion_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.25+0.1*TMath::Log(partE[epart])),scalerange*(0.9));

          h_RH_Reso_pion_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_pion_Erec[ialgo][epart],"L0RMEQ");
          fit_reso_pion_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pion_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[1][ialgo]+2);
          // fit_reso_pion_Erec[ialgo][epart]->Draw("same");
          if(!h_FHCALSTANDALONE_pion_reso_1oErec[icalo][ialgo])h_FHCALSTANDALONE_pion_reso_1oErec[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pion_reso_1oErec_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FHCALSTANDALONE_pion_reso_1oErec[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pion_reso_1oErec[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_pion_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pion_Erec[ialgo][epart]->GetParameter(1));
          if(!h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo])h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pion_meanE_Erec_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_pion_Erec[ialgo][epart]->GetParameter(1));

          // drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          // drawLatexAdd(Form("mean: %1.2f ",fit_reso_pion_Erec[ialgo][epart]->GetParameter(1)),0.15,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          // drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pion_Erec[ialgo][epart]->GetParameter(2)),0.15,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pion_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pion_Erec[ialgo][epart]->GetParameter(1)),0.15,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]){
          float scalerange = 1.0;
          DrawGammaSetMarker(h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[0][ialgo], colorClus[0][ialgo]);
          h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]->Draw("same");
          if(ialgo==kV1)
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 0.15+0.04*TMath::Log(partE[epart]),0.45);
          else if(ialgo==kV3)
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.15+0.1*TMath::Log(partE[epart])),scalerange*(0.85));
          else if(ialgo==k3x3)
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.80));
          else if(ialgo==k5x5)
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC3)
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC5)
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.2+0.07*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kMA)
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.35+0.145*TMath::Log(TMath::Sqrt(partE[epart]))),scalerange*(0.85));
          else
            fit_reso_pionComb_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionComb_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.25+0.1*TMath::Log(partE[epart])),scalerange*(0.9));

          h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_pionComb_Erec[ialgo][epart],"L0RMEQ");
          fit_reso_pionComb_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pionComb_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[0][ialgo]+2);
          fit_reso_pionComb_Erec[ialgo][epart]->Draw("same");
          if(!h_FHCALSTANDALONE_pionComb_reso_1oErec[icalo][ialgo])h_FHCALSTANDALONE_pionComb_reso_1oErec[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pionComb_reso_1oErec_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FHCALSTANDALONE_pionComb_reso_1oErec[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionComb_reso_1oErec[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1));
          if(!h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo])h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pionComb_meanE_Erec_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1)),0.15,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)),0.15,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1)),0.15,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]){
          float scalerange = 1.0;
          DrawGammaSetMarker(h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[2][ialgo], colorClus[2][ialgo]);
          h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]->Draw("same");
          if(ialgo==kV1)
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 3*(0.15+0.04*TMath::Log(partE[epart])),3*(0.45));
          else if(ialgo==kV3)
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 0.8,1.3);
          else if(ialgo==k3x3)
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.80));
          else if(ialgo==k5x5)
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 2*(0.2+0.07*TMath::Log(partE[epart])),2*(0.75));
          else if(ialgo==kC3)
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.1+0.08*TMath::Log(partE[epart])),scalerange*(0.75));
          else if(ialgo==kC5)
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 0.9,1.3);
          else if(ialgo==kMA)
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", 1.4*(0.35+0.135*TMath::Log(TMath::Sqrt(partE[epart]))),1.4*(0.85));
          else
            fit_reso_pionCombCalib_Erec[ialgo][epart] = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s%d_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), ialgo,epart), "gaus", scalerange*(0.25+0.1*TMath::Log(partE[epart])),scalerange*(0.9));

          h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_pionCombCalib_Erec[ialgo][epart],"L0RMEQ");
          fit_reso_pionCombCalib_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pionCombCalib_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[2][ialgo]+2);
          fit_reso_pionCombCalib_Erec[ialgo][epart]->Draw("same");
          if(!h_FHCALSTANDALONE_pionCombCalib_reso_1oErec[icalo][ialgo])h_FHCALSTANDALONE_pionCombCalib_reso_1oErec[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pionCombCalib_reso_1oErec_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FHCALSTANDALONE_pionCombCalib_reso_1oErec[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionCombCalib_reso_1oErec[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1));
          if(!h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo])h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo] 	= new TH1F(Form("h_FHCALSTANDALONE_pionCombCalib_meanE_Erec_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo]->SetBinContent(h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        histEPlotDummy->Draw("same,axis");

        padnum++;
    }
      cPNG->Print(Form("%s/FHCALSTANDALONE_pion_bins_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }


  TCanvas* cReso2x2 = new TCanvas("cReso2x2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2x2, 0.12, 0.01, 0.01, 0.105);
  TH2F* histEPlotDummy_SB_FH                           = new TH2F("histEPlotDummy_SB_FH","histEPlotDummy_SB_FH",1000,0.41, 1.34,1000,0, 0.025);
  SetStyleHistoTH2ForGraphs(histEPlotDummy_SB_FH, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.2);
  histEPlotDummy_SB_FH->GetYaxis()->SetNoExponent();
  histEPlotDummy_SB_FH->GetXaxis()->SetNoExponent();
  histEPlotDummy_SB_FH->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy_SB_FH->GetXaxis()->SetMoreLogLabels(kTRUE);
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
        if(ialgo==kC3 || ialgo==k3x3) continue;
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        if(partE[epart]!=60) continue;
        histEPlotDummy_SB_FH->DrawCopy();
        TLegend* legendEB_wSmear           = GetAndSetLegend2(0.16, 0.95-(2*textSizeLabelsRel), 0.8, 0.95-(5.5*textSizeLabelsRel),1.2*textSizeLabelsPixel, 1, "", 43, 0.15);
        drawLatexAdd(Form("%1.0f GeV < #it{E} < %1.0f GeV",partE[epart],partE[epart+1]),0.16,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("#pi^{-} on %s",str_calorimeter[icalo].Data()),0.93,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart], marker_clusterizer[0], 3, colorClus[0][ialgo], colorClus[0][ialgo]);
          h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_pionComb_Erec[ialgo][epart]    = new TF1(Form("fit_reso_pionComb_Erec_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                                      0.41, 0.9);
          fit_reso_pionComb_Erec[ialgo][epart]->SetParameter(0, 0.015);
          fit_reso_pionComb_Erec[ialgo][epart]->SetParameter(1, 0.73);
          fit_reso_pionComb_Erec[ialgo][epart]->SetParameter(2, 0.08);
          fit_reso_pionComb_Erec[ialgo][epart]->SetParameter(3, 0.03);
          h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_pionComb_Erec[ialgo][epart],"L0MEQ", "", 0.41, 0.9);
          fit_reso_pionComb_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pionComb_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 3, colorClus[0][ialgo]-1);
          fit_reso_pionComb_Erec[ialgo][epart]->Draw("same");

          legendEB_wSmear->AddEntry(h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart], Form("uncalibrated: mean: %1.1f, #sigma/#it{E}: %1.2f",fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1)), "pl");
          // drawLatexAdd("uncalibrated:",0.16,0.84,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          // drawLatexAdd(Form("mean: %1.2f ",fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1)),0.16,0.78,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          // drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)),0.16,0.72,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1)),0.16,0.66,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart], marker_clusterizer[2], 3, colorClus[2][ialgo], colorClus[2][ialgo]);
          h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_pionCombCalib_Erec[ialgo][epart]    = new TF1(Form("fit_reso_pionCombCalib_Erec_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                                      0.7, 1.3);
          fit_reso_pionCombCalib_Erec[ialgo][epart]->SetParameter(0, 0.015);
          fit_reso_pionCombCalib_Erec[ialgo][epart]->SetParameter(1, 1.0);
          fit_reso_pionCombCalib_Erec[ialgo][epart]->SetParameter(2, 0.11);
          fit_reso_pionCombCalib_Erec[ialgo][epart]->SetParameter(3, 0.05);
          h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_pionCombCalib_Erec[ialgo][epart],"L0MEQ", "", 0.7, 1.3);
          fit_reso_pionCombCalib_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pionCombCalib_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 3, colorClus[2][ialgo]+1);
          fit_reso_pionCombCalib_Erec[ialgo][epart]->Draw("same");
          
          legendEB_wSmear->AddEntry(h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart], Form("calibrated: mean: %1.1f, #sigma/#it{E}: %1.2f",fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1)), "pl");

          // drawLatexAdd("calibrated:",0.93,0.84,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          // drawLatexAdd(Form("mean: %1.2f ",fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1)),0.93,0.78,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          // drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(2)),0.93,0.72,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1)),0.93,0.66,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        if(h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart], marker_clusterizer[1], 3, colorClus[1][ialgo], colorClus[1][ialgo]);
          h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]    = new TF1(Form("fit_reso_pionCombCalib_smeared_Erec_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                                      0.7, 1.3);
          fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->SetParameter(0, 0.015);
          fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->SetParameter(1, 1.0);
          fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->SetParameter(2, 0.11);
          fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->SetParameter(3, 0.05);
          h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart]->Fit(fit_reso_pionCombCalib_smeared_Erec[ialgo][epart],"L0MEQ", "", 0.7, 1.3);
          fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_pionCombCalib_smeared_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 3, colorClus[1][ialgo]+1);
          fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->Draw("same");

          legendEB_wSmear->AddEntry(h_RH_Reso_pionCombCalib_bin_smeared_Erec[icalo][ialgo][epart], Form("calib+smear: mean: %1.1f, #sigma/#it{E}: %1.2f",fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->GetParameter(1)), "pl");

          // drawLatexAdd("calibrated:",0.93,0.84,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          // drawLatexAdd(Form("mean: %1.2f ",fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->GetParameter(1)),0.93,0.78,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          // drawLatexAdd(Form("sigma: %1.2f ",fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->GetParameter(2)),0.93,0.72,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib_smeared_Erec[ialgo][epart]->GetParameter(1)),0.93,0.66,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
    legendEB_wSmear->Draw();
    }
    histEPlotDummy_SB_FH->Draw("same,axis");
    cReso2x2->Print(Form("%s/FHCALSTANDALONE_pion_SingleBin_smear_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
  
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        if(partE[epart]!=60) continue;
        histEPlotDummy_SB_FH->DrawCopy();
        TLegend* legendEB_wSmear           = GetAndSetLegend2(0.16, 0.95-(2*textSizeLabelsRel), 0.8, 0.95-(5.5*textSizeLabelsRel),1.2*textSizeLabelsPixel, 1, "", 43, 0.15);
        drawLatexAdd(Form("%1.0f GeV < #it{E} < %1.0f GeV",partE[epart],partE[epart+1]),0.16,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("#pi^{-} on %s",str_calorimeter[icalo].Data()),0.93,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]){
          h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_pionComb_Erec[ialgo][epart]->Draw("same");

          legendEB_wSmear->AddEntry(h_RH_Reso_pionComb_bin_Erec[icalo][ialgo][epart], Form("uncalibrated: mean: %1.1f, #sigma/#it{E}: %1.2f",fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1)), "pl");
        }
        if(h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]){
          h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_pionCombCalib_Erec[ialgo][epart]->Draw("same");
          
          legendEB_wSmear->AddEntry(h_RH_Reso_pionCombCalib_bin_Erec[icalo][ialgo][epart], Form("calibrated: mean: %1.1f, #sigma/#it{E}: %1.2f",fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionCombCalib_Erec[ialgo][epart]->GetParameter(1)), "pl");
        }
          legendEB_wSmear->AddEntry((TObject*)0, " ", "");
    legendEB_wSmear->Draw();
    }
    histEPlotDummy_SB_FH->Draw("same,axis");
      cReso2x2->Print(Form("%s/FHCALSTANDALONE_pion_SingleBin_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }

  // cReso2x2 = new TCanvas("cReso2x2","",0,0,1100,1000);
  // DrawGammaCanvasSettings( cReso2x2, 0.12, 0.01, 0.01, 0.105);
  // TH2F* histEPlotDummy_SingBin                           = new TH2F("histEPlotDummy_SingBin","histEPlotDummy_SingBin",1000,0.61, 1.14,1000,0, 0.21);
  // SetStyleHistoTH2ForGraphs(histEPlotDummy_SingBin, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.2);
  // histEPlotDummy_SingBin->GetXaxis()->SetNoExponent();
  // histEPlotDummy_SingBin->GetYaxis()->SetNdivisions(505,kTRUE);
  // histEPlotDummy_SingBin->GetXaxis()->SetMoreLogLabels(kTRUE);
  // for(int icalo=0;icalo<_active_calo;icalo++){
  //   if(icalo!=1)continue;
  //   for(int ialgo=0;ialgo<_active_algo;ialgo++){
  //       if(ialgo==kC3 || ialgo==k3x3) continue;
  //     for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
  //       if(partE[epart]!=40) continue;
  //       histEPlotDummy_SingBin->DrawCopy();
  //       drawLatexAdd(Form("%1.0f GeV < #it{E} < %1.0f GeV",partE[epart],partE[epart+1]),0.16,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  //       drawLatexAdd(Form("#gamma on %s",str_calorimeter[icalo].Data()),0.93,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  //       if(h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]){
  //         float scalerange = 1;
  //         h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]->Draw("same");
  //         fit_reso_gamma_Erec[ialgo][epart]->Draw("same");

  //         drawLatexAdd("uncalibrated:",0.16,0.84,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  //         drawLatexAdd(Form("mean: %1.2f ",fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1)),0.16,0.78,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  //         drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)),0.16,0.72,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  //         drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1)),0.16,0.66,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
  //       }
  //       if(h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]){
  //         h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]->Draw("same");
  //         fit_reso_gammaCalib_Erec[ialgo][epart]->Draw("same");

  //         drawLatexAdd("calibrated:",0.93,0.84,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  //         drawLatexAdd(Form("mean: %1.2f ",fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1)),0.93,0.78,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  //         drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(2)),0.93,0.72,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  //         drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1)),0.93,0.66,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  //       }
  //   }
  //     histEPlotDummy_SingBin->Draw("same,axis");
  //     cReso2x2->Print(Form("%s/FEMCSTANDALONE_gamma_SingleBin_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
  //   }
  // }

	split_canvas(cPNG, "cPNG1", nEne);

  histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0.6, 1.2,1000,0, 0.11);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  TF1 *fit_reso_gamma_Etrue[nEne] = {NULL};
  TF1 *fit_reso_gammaCalib_Etrue[nEne] = {NULL};
  // TH1F*  h_gamma_resolution_femc_E 	= (TH1F*)h_RH_Reso_bin_E[icalo][ialgo][epart]->Clone("h_gamma_resolution_femc_E");
  // for(int i=1;i<h_gamma_resolution_femc_E->GetNbinsX();i++) h_gamma_resolution_femc_E->SetBinContent(i,0);
  TH1F*  h_FEMCSTANDALONE_gamma_reso_1oE[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FEMCSTANDALONE_gamma_meanE_E[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FEMCSTANDALONE_gammaCalib_reso_1oE[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FEMCSTANDALONE_gammaCalib_meanE_E[_active_calo][_active_algo] 	= {NULL};
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=1)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      Int_t padnum =0;
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1);

        histEPlotDummy->DrawCopy();
        if(h_RH_Reso_gamma_bin_E[icalo][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_gamma_bin_E[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[0][ialgo], colorClus[0][ialgo]);
          h_RH_Reso_gamma_bin_E[icalo][ialgo][epart]->Draw("same");
          fit_reso_gamma_Etrue[epart]    = new TF1(Form("fit_reso_gamma_Etrue_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                                      0.60, 0.86);
          fit_reso_gamma_Etrue[epart]->SetParameter(0, 0.06);
          fit_reso_gamma_Etrue[epart]->SetParameter(1, 0.82);
          fit_reso_gamma_Etrue[epart]->SetParameter(2, 0.01);
          fit_reso_gamma_Etrue[epart]->SetParameter(3, 0.012);
          if(ialgo==kV3){
            fit_reso_gamma_Etrue[epart]->SetParameter(0, 0.03);
            fit_reso_gamma_Etrue[epart]->SetParameter(1, 0.77);
            h_RH_Reso_gamma_bin_E[icalo][ialgo][epart]->Fit(fit_reso_gamma_Etrue[epart],"L0MEQ", "", 0.60, 0.85);
          } else {
            h_RH_Reso_gamma_bin_E[icalo][ialgo][epart]->Fit(fit_reso_gamma_Etrue[epart],"L0MEQ", "", 0.60, 0.88);
          }
          fit_reso_gamma_Etrue[epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_gamma_Etrue[epart], linestyle_clusterizer[ialgo], 2, colorClus[0][ialgo]+2);
          fit_reso_gamma_Etrue[epart]->Draw("same");
          if(!h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo])h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo] 	= new TH1F(Form("h_gamma_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_gamma_Etrue[epart]->GetParameter(2)/fit_reso_gamma_Etrue[epart]->GetParameter(1));
          if(!h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo])h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo] 	= new TH1F(Form("h_FEMCSTANDALONE_gamma_meanE_E_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_gamma_Etrue[epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_gamma_Etrue[epart]->GetParameter(1)),0.15,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gamma_Etrue[epart]->GetParameter(2)),0.15,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gamma_Etrue[epart]->GetParameter(2)/fit_reso_gamma_Etrue[epart]->GetParameter(1)),0.15,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_gammaCalib_bin_E[icalo][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_gammaCalib_bin_E[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[2][ialgo], colorClus[2][ialgo]);
          h_RH_Reso_gammaCalib_bin_E[icalo][ialgo][epart]->Draw("same");

          fit_reso_gammaCalib_Etrue[epart]    = new TF1(Form("fit_reso_gammaCalib_Etrue_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                            0.9, 1.1);
          fit_reso_gammaCalib_Etrue[epart]->SetParameter(0, 0.06);
          fit_reso_gammaCalib_Etrue[epart]->SetParameter(1, 1.0);
          fit_reso_gammaCalib_Etrue[epart]->SetParameter(2, 0.01);
          fit_reso_gammaCalib_Etrue[epart]->SetParameter(3, 0.012);
          h_RH_Reso_gammaCalib_bin_E[icalo][ialgo][epart]->Fit(fit_reso_gammaCalib_Etrue[epart],"L0MEQ", "", 0.8, 1.1);

          fit_reso_gammaCalib_Etrue[epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_gammaCalib_Etrue[epart], linestyle_clusterizer[ialgo], 2, colorClus[2][ialgo]+2);
          fit_reso_gammaCalib_Etrue[epart]->Draw("same");
          if(!h_FEMCSTANDALONE_gammaCalib_reso_1oE[icalo][ialgo])h_FEMCSTANDALONE_gammaCalib_reso_1oE[icalo][ialgo] 	= new TH1F(Form("h_gammaCalib_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FEMCSTANDALONE_gammaCalib_reso_1oE[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gammaCalib_reso_1oE[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_gammaCalib_Etrue[epart]->GetParameter(2)/fit_reso_gammaCalib_Etrue[epart]->GetParameter(1));
          if(!h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo])h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo] 	= new TH1F(Form("h_FEMCSTANDALONE_gammaCalib_meanE_E_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_gammaCalib_Etrue[epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_gammaCalib_Etrue[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gammaCalib_Etrue[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gammaCalib_Etrue[epart]->GetParameter(2)/fit_reso_gammaCalib_Etrue[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
      histEPlotDummy->Draw("same,axis");
      padnum++;
    }
      // cPNG->Print(Form("%s/FEMCSTANDALONE_gamma_bins_%s%s_Etrue.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }

  TF1 *fit_reso_gamma_Erec[_active_algo][nEne] = {NULL};
  TH1F*  h_FEMCSTANDALONE_gamma_reso_1oErec[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FEMCSTANDALONE_gamma_meanE_Erec[_active_calo][_active_algo] 	= {NULL};
  TF1 *fit_reso_gammaCalib_Erec[_active_algo][nEne] = {NULL};
  TH1F*  h_FEMCSTANDALONE_gammaCalib_reso_1oErec[_active_calo][_active_algo] 	= {NULL};
  TH1F*  h_FEMCSTANDALONE_gammaCalib_meanE_Erec[_active_calo][_active_algo] 	= {NULL};
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=1)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      Int_t padnum =0;
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        DrawVirtualPadSettings( cPNG->cd(padnum+1), 0.1, 0.02, 0.015, 0.115);
        cPNG->cd(padnum+1);

        histEPlotDummy->DrawCopy();
        if(h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[0][ialgo], colorClus[0][ialgo]);
          h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_gamma_Erec[ialgo][epart]    = new TF1(Form("fit_reso_gamma_Erec_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                                      0.60, 0.86);
          fit_reso_gamma_Erec[ialgo][epart]->SetParameter(0, 0.06);
          fit_reso_gamma_Erec[ialgo][epart]->SetParameter(1, 0.82);
          fit_reso_gamma_Erec[ialgo][epart]->SetParameter(2, 0.01);
          fit_reso_gamma_Erec[ialgo][epart]->SetParameter(3, 0.012);
          if(ialgo==kV3){
            fit_reso_gamma_Erec[ialgo][epart]->SetParameter(0, 0.03);
            fit_reso_gamma_Erec[ialgo][epart]->SetParameter(1, 0.77);
            h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_gamma_Erec[ialgo][epart],"L0MEQ", "", 0.60, 0.85);
          } else {
            h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_gamma_Erec[ialgo][epart],"L0MEQ", "", 0.60, 0.88);
          }
          fit_reso_gamma_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_gamma_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[0][ialgo]+2);
          fit_reso_gamma_Erec[ialgo][epart]->Draw("same");
          if(!h_FEMCSTANDALONE_gamma_reso_1oErec[icalo][ialgo])h_FEMCSTANDALONE_gamma_reso_1oErec[icalo][ialgo] 	= new TH1F(Form("h_gamma_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FEMCSTANDALONE_gamma_reso_1oErec[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gamma_reso_1oErec[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1));
          if(!h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo])h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo] 	= new TH1F(Form("h_FEMCSTANDALONE_gamma_meanE_E_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1)),0.15,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)),0.15,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1)),0.15,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart], marker_clusterizer[ialgo], 3, colorClus[2][ialgo], colorClus[0][ialgo]);
          h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_gammaCalib_Erec[ialgo][epart]    = new TF1(Form("fit_reso_gammaCalib_Erec_%s%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),epart),"(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                                      0.9, 1.1);
          fit_reso_gammaCalib_Erec[ialgo][epart]->SetParameter(0, 0.06);
          fit_reso_gammaCalib_Erec[ialgo][epart]->SetParameter(1, 1.0);
          fit_reso_gammaCalib_Erec[ialgo][epart]->SetParameter(2, 0.01);
          fit_reso_gammaCalib_Erec[ialgo][epart]->SetParameter(3, 0.012);
          h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]->Fit(fit_reso_gammaCalib_Erec[ialgo][epart],"L0MEQ", "", 0.8, 1.1);

          fit_reso_gammaCalib_Erec[ialgo][epart]->SetNpx(10000);
          DrawGammaSetMarkerTF1(fit_reso_gammaCalib_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorClus[2][ialgo]+2);
          fit_reso_gammaCalib_Erec[ialgo][epart]->Draw("same");
          if(!h_FEMCSTANDALONE_gammaCalib_reso_1oErec[icalo][ialgo])h_FEMCSTANDALONE_gammaCalib_reso_1oErec[icalo][ialgo] 	= new TH1F(Form("h_gammaCalib_resolution_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 1.0);
          h_FEMCSTANDALONE_gammaCalib_reso_1oErec[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gammaCalib_reso_1oErec[icalo][ialgo]->FindBin(1./TMath::Sqrt((partE[epart]+partE[epart+1])/2)),100*fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1));
          if(!h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo])h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo] 	= new TH1F(Form("h_FEMCSTANDALONE_gammaCalib_meanE_E_%s%s_1oE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500, 0.0, 100.0);
          h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo]->SetBinContent(h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo]->FindBin((partE[epart]+partE[epart+1])/2),fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1));

          drawLatexAdd(Form("%1.1f GeV < #it{E} < %1.1f GeV #pi^{-} on %s",partE[epart],partE[epart+1],str_calorimeter[icalo].Data()),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

          drawLatexAdd(Form("mean: %1.2f ",fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        padnum++;
    }
      cPNG->Print(Form("%s/FEMCSTANDALONE_gamma_bins_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }
  cReso2x2 = new TCanvas("cReso2x2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2x2, 0.12, 0.01, 0.01, 0.105);
  TH2F* histEPlotDummy_SingBin                           = new TH2F("histEPlotDummy_SingBin","histEPlotDummy_SingBin",1000,0.61, 1.14,1000,0, 0.139);
  SetStyleHistoTH2ForGraphs(histEPlotDummy_SingBin, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.2);
  histEPlotDummy_SingBin->GetXaxis()->SetNoExponent();
  histEPlotDummy_SingBin->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy_SingBin->GetXaxis()->SetMoreLogLabels(kTRUE);
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=1)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
        if(ialgo==kC3 || ialgo==k3x3) continue;
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        if(partE[epart]!=40) continue;
        histEPlotDummy_SingBin->DrawCopy();
        drawLatexAdd(Form("%1.0f GeV < #it{E} < %1.0f GeV",partE[epart],partE[epart+1]),0.16,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("#gamma on %s",str_calorimeter[icalo].Data()),0.93,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]){
          float scalerange = 1;
          h_RH_Reso_gamma_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_gamma_Erec[ialgo][epart]->Draw("same");

          drawLatexAdd("uncalibrated:",0.16,0.84,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("mean: %1.2f ",fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1)),0.16,0.78,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)),0.16,0.72,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1)),0.16,0.66,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        }
        if(h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]){
          h_RH_Reso_gammaCalib_bin_Erec[icalo][ialgo][epart]->Draw("same");
          fit_reso_gammaCalib_Erec[ialgo][epart]->Draw("same");

          drawLatexAdd("calibrated:",0.93,0.84,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("mean: %1.2f ",fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1)),0.93,0.78,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("sigma: %1.2f ",fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(2)),0.93,0.72,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
          drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gammaCalib_Erec[ialgo][epart]->GetParameter(1)),0.93,0.66,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
    }
      histEPlotDummy_SingBin->Draw("same,axis");
      cReso2x2->Print(Form("%s/FEMCSTANDALONE_gamma_SingleBin_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }

  TH2F* histEPlotDummy_SingBin2                           = new TH2F("histEPlotDummy_SingBin2","histEPlotDummy_SingBin2",1000,0.21, 1.04,1000,0, 0.119);
  SetStyleHistoTH2ForGraphs(histEPlotDummy_SingBin2, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,1.2);
  histEPlotDummy_SingBin2->GetYaxis()->SetNoExponent();
  histEPlotDummy_SingBin2->GetXaxis()->SetNoExponent();
  histEPlotDummy_SingBin2->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy_SingBin2->GetXaxis()->SetMoreLogLabels(kTRUE);
  if(1){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
    TLegend* legendHE_pion           = GetAndSetLegend2(0.16, 0.95-(3*textSizeLabelsRel), 0.75, 0.95-(6*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);
      if(ialgo==kC3 || ialgo==k3x3) continue;
      histEPlotDummy_SingBin2->DrawCopy();
      for(Int_t epart=0; epart<NELEMS(partE)-1; epart++){
        if(partE[epart]!=60) continue;
        drawLatexAdd(Form("%s clusterizer: E_{seed} = 500 MeV, E_{agg} = 100 MeV",str_clusterizer[ialgo].Data()),0.16,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        drawLatexAdd(Form("%1.0f GeV < #it{E}_{part} < %1.0f GeV",partE[epart],partE[epart+1]),0.16,0.85,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
        // drawLatexAdd(Form("#pi^{-} on %s",str_calorimeter[0].Data()),0.93,0.90,1.1*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(h_RH_Reso_pionComb_bin_Erec[0][ialgo][epart]){
          h_RH_Reso_pionComb_bin_Erec[0][ialgo][epart]->Scale(2.5);
          DrawGammaSetMarker(h_RH_Reso_pionComb_bin_Erec[0][ialgo][epart], markerStylePID[2]-4, 1.5*markerSizePID[2], colorPID[2], colorPID[2]);
          h_RH_Reso_pionComb_bin_Erec[0][ialgo][epart]->Draw("same");

          DrawGammaSetMarkerTF1(fit_reso_pionComb_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorPID[2]+2);
          legendHE_pion->AddEntry(h_RH_Reso_pionComb_bin_Erec[0][ialgo][epart],Form("#pi^{-} in FHCAL+FEMC: mean: %1.2f, #sigma/#it{E}: %1.2f %%",fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pionComb_Erec[ialgo][epart]->GetParameter(1)),"pl");
          fit_reso_pionComb_Erec[ialgo][epart]=ScaleTF1(fit_reso_pionComb_Erec[ialgo][epart],2.5,Form("scaledfunc%d%d",ialgo,epart));
          DrawGammaSetMarkerTF1(fit_reso_pionComb_Erec[ialgo][epart], linestyle_clusterizer[ialgo], 2, colorPID[2]+2);
          fit_reso_pionComb_Erec[ialgo][epart]->Draw("same");

        }
        if(h_RH_Reso_pion_bin_Erec[0][ialgo][epart]){
          h_RH_Reso_pion_bin_Erec[0][ialgo][epart]->Scale(2.5);
          DrawGammaSetMarker(h_RH_Reso_pion_bin_Erec[0][ialgo][epart], markerStylePID[1]-4, 1.5*markerSizePID[1], colorPID[1], colorPID[1]);
          h_RH_Reso_pion_bin_Erec[0][ialgo][epart]->Draw("same");

          DrawGammaSetMarkerTF1(fit_reso_pion_Erec[ialgo][epart], linestyle_clusterizer[ialgo]+1, 2, colorPID[1]+2);
          legendHE_pion->AddEntry(h_RH_Reso_pion_bin_Erec[0][ialgo][epart],Form("#pi^{-} in FHCAL: mean: %1.2f, #sigma/#it{E}: %1.2f %%",fit_reso_pion_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_pion_Erec[ialgo][epart]->GetParameter(2)/fit_reso_pion_Erec[ialgo][epart]->GetParameter(1)),"pl");
          fit_reso_pion_Erec[ialgo][epart]=ScaleTF1(fit_reso_pion_Erec[ialgo][epart],2.5,Form("scaledfunc%d%d",ialgo,epart));
          DrawGammaSetMarkerTF1(fit_reso_pion_Erec[ialgo][epart], linestyle_clusterizer[ialgo]+1, 2, colorPID[1]+2);
          fit_reso_pion_Erec[ialgo][epart]->Draw("same");

        }
        if(h_RH_Reso_gamma_bin_Erec[1][ialgo][epart]){
          DrawGammaSetMarker(h_RH_Reso_gamma_bin_Erec[1][ialgo][epart], markerStylePID[3]-4, 1.5*markerSizePID[3], colorPID[3], colorPID[3]);
          h_RH_Reso_gamma_bin_Erec[1][ialgo][epart]->Draw("same");
          DrawGammaSetMarkerTF1(fit_reso_gamma_Erec[ialgo][epart], linestyle_clusterizer[ialgo]-1, 2, colorPID[3]+2);
          fit_reso_gamma_Erec[ialgo][epart]->Draw("same");

          legendHE_pion->AddEntry(h_RH_Reso_gamma_bin_Erec[1][ialgo][epart],Form("#gamma in FEMC: mean: %1.2f, #sigma/#it{E}: %1.2f %%",fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1),100*fit_reso_gamma_Erec[ialgo][epart]->GetParameter(2)/fit_reso_gamma_Erec[ialgo][epart]->GetParameter(1)),"pl");
        }
    }
    legendHE_pion->Draw();
      histEPlotDummy_SingBin2->Draw("same,axis");
      cReso2x2->Print(Form("%s/CALO_ExampleSingleBin_%s_Erec.%s", outputDir.Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }
	split_canvas(cPNG, "cPNG1", nEne);

// RESOLUTION PLOT 2

  TCanvas* cReso2 = new TCanvas("cReso2","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2, 0.095, 0.01, 0.01, 0.105);
  TH2F* histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,0, 0.49,1000,0, 54);
  SetStyleHistoTH2ForGraphs(histResoDummy, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      histResoDummy->Draw();
        DrawGammaSetMarker(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo], 46, 3, kBlack, kBlack);

      h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->Draw("same,p");
      TF1* fit_reso_pion_1oE = new TF1("fit_reso_pion_1oE", "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
      fit_reso_pion_1oE->SetParameter(2,0.3);
      fit_reso_pion_1oE->FixParameter(2,0.08);
      h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->Fit(fit_reso_pion_1oE,"QRMNE");
      fit_reso_pion_1oE->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_pion_1oE, 3, 3, kRed+2);
      fit_reso_pion_1oE->Draw("same");
      drawLatexAdd(Form("#pi^{-} in %s",str_calorimeter[icalo].Data()),0.90,0.25,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd(Form("#frac{#sigma}{#it{E}} =  #frac{%1.2f}{#sqrt{#it{E}}} #oplus %1.2f",fit_reso_pion_1oE->GetParameter(1),fit_reso_pion_1oE->GetParameter(0)),0.90,0.35,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


      // cReso2->Print(Form("%s/FHCALSTANDALONE_pion_reso_1oE_%s%s_true.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
    }
  }

  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
    histResoDummy->Draw();
    TLegend* legendHE_pion           = GetAndSetLegend2(0.1, 0.95, 0.7, 0.95-(7*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

    TF1* fit_reso_clusdep_pion_1oE[_active_algo];
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      DrawGammaSetMarker(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo], marker_clusterizer[ialgo], 3, colorClus[0][ialgo], colorClus[0][ialgo]);
      h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->SetLineColor(colorClus[0][ialgo]-10);
      h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
      h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->SetLineWidth(3);
      h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->Draw("same,p");
      fit_reso_clusdep_pion_1oE[ialgo] = new TF1(Form("fit_reso_clusdep_pion_1oE_%d%d",ialgo,icalo), "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
      fit_reso_clusdep_pion_1oE[ialgo]->SetParameter(2,0.3);
      fit_reso_clusdep_pion_1oE[ialgo]->FixParameter(2,0.08);
      h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->Fit(fit_reso_clusdep_pion_1oE[ialgo],"QRMNE");
      // fit_reso_clusdep_pion_1oE[ialgo]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_clusdep_pion_1oE[ialgo], 3, 3, colorClus[0][ialgo]-10);
      fit_reso_clusdep_pion_1oE[ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
      fit_reso_clusdep_pion_1oE[ialgo]->Draw("same");
      legendHE_pion->AddEntry(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo],Form("%s: #sigma/#it{E} =  %1.2f/#sqrt{#it{E}} #oplus %1.2f",str_clusterizer[ialgo].Data(),fit_reso_clusdep_pion_1oE[ialgo]->GetParameter(1),fit_reso_clusdep_pion_1oE[ialgo]->GetParameter(0)),"pl");
    }
    legendHE_pion->Draw();
    // drawLatexAdd(Form("#pi^{-} in %s",str_calorimeter[icalo].Data()),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#pi^{-} in %s",str_calorimeter[icalo].Data()),0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("placeholder data"),0.5,0.5,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE,kOrange+2);


    cReso2->Print(Form("%s/FHCALSTANDALONE_pion_reso_1oE_%s_diffCLSRZR_true.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
  }

    if(1){
      histResoDummy->Draw();
      TLegend* legendHE_pion           = GetAndSetLegend2(0.15, 0.95-(2*textSizeLabelsRel), 0.75, 0.95-(6*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

      TF1* fit_reso_clusdep_pion_1oE[_active_algo];
      TF1* fit_reso_clusdep_pionComb_1oE[_active_algo];
      for(int icalo=0;icalo<1;icalo++){
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          if(!(ialgo==kMA)) continue;
          if(h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]){
            DrawGammaSetMarker(h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo], markerStylePID[2]-4, 1.5*markerSizePID[2], colorPID[2], colorPID[2]);
            h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]->SetLineColor(colorPID[2]-10);
            h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
            h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]->SetLineWidth(3);
            h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]->Draw("same,p");
            fit_reso_clusdep_pionComb_1oE[ialgo] = new TF1(Form("fit_reso_clusdep_pionComb_1oE_%d%d",ialgo,icalo), "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
            fit_reso_clusdep_pionComb_1oE[ialgo]->SetParameter(2,0.3);
            fit_reso_clusdep_pionComb_1oE[ialgo]->FixParameter(2,0.08);
            h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo]->Fit(fit_reso_clusdep_pionComb_1oE[ialgo],"QRMNE");
            // fit_reso_clusdep_pionComb_1oE[ialgo]->SetNpx(10000);
            DrawGammaSetMarkerTF1(fit_reso_clusdep_pionComb_1oE[ialgo], 3, 3, colorPID[2]-10);
            fit_reso_clusdep_pionComb_1oE[ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
            fit_reso_clusdep_pionComb_1oE[ialgo]->Draw("same");
            legendHE_pion->AddEntry(h_FHCALSTANDALONE_pionComb_reso_1oE[icalo][ialgo],Form("#pi^{-} in FHCAL (with material): #sigma/#it{E} =  %1.2f/#sqrt{#it{E}} #oplus %1.2f",fit_reso_clusdep_pionComb_1oE[ialgo]->GetParameter(1),fit_reso_clusdep_pionComb_1oE[ialgo]->GetParameter(0)),"pl");
          }
        }
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          if(!(ialgo==kMA)) continue;
          if(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]){
            DrawGammaSetMarker(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo], markerStylePID[1]-4, 1.5*markerSizePID[1], colorPID[1], colorPID[1]);
            h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->SetLineColor(colorPID[1]-10);
            h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
            h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->SetLineWidth(3);
            h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->Draw("same,p");
            fit_reso_clusdep_pion_1oE[ialgo] = new TF1(Form("fit_reso_clusdep_pion_1oE_%d%d",ialgo,icalo), "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
            fit_reso_clusdep_pion_1oE[ialgo]->SetParameter(2,0.3);
            fit_reso_clusdep_pion_1oE[ialgo]->FixParameter(2,0.08);
            h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo]->Fit(fit_reso_clusdep_pion_1oE[ialgo],"QRMNE");
            // fit_reso_clusdep_pion_1oE[ialgo]->SetNpx(10000);
            DrawGammaSetMarkerTF1(fit_reso_clusdep_pion_1oE[ialgo], 3, 3, colorPID[1]-10);
            fit_reso_clusdep_pion_1oE[ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
            fit_reso_clusdep_pion_1oE[ialgo]->Draw("same");
            // legendHE_pion->AddEntry(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo],Form("%s: #sigma/#it{E} =  %1.2f/#sqrt{#it{E}} #oplus %1.2f",str_clusterizer[ialgo].Data(),fit_reso_clusdep_pion_1oE[ialgo]->GetParameter(1),fit_reso_clusdep_pion_1oE[ialgo]->GetParameter(0)),"pl");
            legendHE_pion->AddEntry(h_FHCALSTANDALONE_pion_reso_1oE[icalo][ialgo],Form("#pi^{-} in FHCAL (w/o material): #sigma/#it{E} =  %1.2f/#sqrt{#it{E}} #oplus %1.2f",fit_reso_clusdep_pion_1oE[ialgo]->GetParameter(1),fit_reso_clusdep_pion_1oE[ialgo]->GetParameter(0)),"pl");
          }
        }
      }
      for(int icalo=1;icalo<2;icalo++){
        TF1* fit_reso_clusdep_gamma_1oE[_active_algo];
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          if(!(ialgo==kMA)) continue;
          if(h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]){
            DrawGammaSetMarker(h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo], markerStylePID[3]-4, 1.5*markerSizePID[3], colorPID[3], colorPID[3]);
            h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->SetLineColor(colorPID[3]-10);
            h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
            h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->SetLineWidth(3);
            h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->Draw("same,p");
            fit_reso_clusdep_gamma_1oE[ialgo] = new TF1(Form("fit_reso_clusdep_gamma_1oE_%d%d",ialgo,icalo), "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
            fit_reso_clusdep_gamma_1oE[ialgo]->SetParameter(2,0.3);
            fit_reso_clusdep_gamma_1oE[ialgo]->FixParameter(2,0.08);
            h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->Fit(fit_reso_clusdep_gamma_1oE[ialgo],"QRMNE");
            // fit_reso_clusdep_gamma_1oE[ialgo]->SetNpx(10000);
            DrawGammaSetMarkerTF1(fit_reso_clusdep_gamma_1oE[ialgo], 3, 3, colorPID[3]-10);
            fit_reso_clusdep_gamma_1oE[ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
            fit_reso_clusdep_gamma_1oE[ialgo]->Draw("same");
            legendHE_pion->AddEntry(h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo],Form("#gamma in FEMC: #sigma/#it{E} =  %1.2f/#sqrt{#it{E}} #oplus %1.2f",fit_reso_clusdep_gamma_1oE[ialgo]->GetParameter(1),fit_reso_clusdep_gamma_1oE[ialgo]->GetParameter(0)),"pl");
          }
        }
      }
      legendHE_pion->Draw();
      drawLatexAdd("MA clusterizer: E_{seed} = 500 MeV, E_{agg} = 100 MeV",0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
      // drawLatexAdd(Form("#pi^{-} in %s",str_calorimeter[icalo].Data()),0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      // drawLatexAdd(Form("placeholder data"),0.5,0.5,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE,kOrange+2);
      cReso2->Print(Form("%s/CALOCOMPARISON_reso_1oE_true.%s", outputDir.Data(), suffix.Data()));
    }
  TH2F* histResoDummyGamma                           = new TH2F("histResoDummyGamma","histResoDummyGamma",1000,0, 0.49,1000,0, 11.9);
  SetStyleHistoTH2ForGraphs(histResoDummyGamma, "1 / #sqrt{#it{E} (GeV)}","#sigma / #it{E} (%)", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histResoDummyGamma->GetXaxis()->SetNoExponent();
  histResoDummyGamma->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummyGamma->GetXaxis()->SetMoreLogLabels(kTRUE);
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=1)continue;
    histResoDummyGamma->Draw();
    TLegend* legendHE_gamma           = GetAndSetLegend2(0.1, 0.95, 0.7, 0.95-(7*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);

    TF1* fit_reso_clusdep_gamma_1oE[_active_algo];
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      DrawGammaSetMarker(h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo], marker_clusterizer[ialgo], 3, colorClus[0][ialgo], colorClus[0][ialgo]);
      h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->SetLineColor(colorClus[0][ialgo]-10);
      h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
      h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->SetLineWidth(3);
      h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->Draw("same,p");
      fit_reso_clusdep_gamma_1oE[ialgo] = new TF1(Form("fit_reso_clusdep_gamma_1oE_%d%d",ialgo,icalo), "[0]+[1]*x+[2]*TMath::Power(x,2)", 0,0.5);
      fit_reso_clusdep_gamma_1oE[ialgo]->SetParameter(2,0.3);
      fit_reso_clusdep_gamma_1oE[ialgo]->FixParameter(2,0.08);
      h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo]->Fit(fit_reso_clusdep_gamma_1oE[ialgo],"QRMNE");
      // fit_reso_clusdep_gamma_1oE[ialgo]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso_clusdep_gamma_1oE[ialgo], 3, 3, colorClus[0][ialgo]-10);
      fit_reso_clusdep_gamma_1oE[ialgo]->SetLineStyle(linestyle_clusterizer[ialgo]);
      fit_reso_clusdep_gamma_1oE[ialgo]->Draw("same");
      legendHE_gamma->AddEntry(h_FEMCSTANDALONE_gamma_reso_1oE[icalo][ialgo],Form("%s: #sigma/#it{E} =  %1.2f/#sqrt{#it{E}} #oplus %1.2f",str_clusterizer[ialgo].Data(),fit_reso_clusdep_gamma_1oE[ialgo]->GetParameter(1),fit_reso_clusdep_gamma_1oE[ialgo]->GetParameter(0)),"pl");
    }
    legendHE_gamma->Draw();
    // drawLatexAdd(Form("#pi^{-} in %s",str_calorimeter[icalo].Data()),0.15,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("#gamma in %s",str_calorimeter[icalo].Data()),0.95,0.15,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("placeholder data"),0.5,0.5,1.1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE,kOrange+2);


    cReso2->Print(Form("%s/FEMCSTANDALONE_gamma_reso_1oE_%s_diffCLSRZR_true.%s", outputDir.Data(),str_calorimeter[icalo].Data(), suffix.Data()));
  }


  TCanvas* cReso2x = new TCanvas("cReso2x","",0,0,1100,1000);
  DrawGammaCanvasSettings( cReso2x, 0.1, 0.01, 0.01, 0.105);
  cReso2x->SetLogx();
  histResoDummy                           = new TH2F("histResoDummy","histResoDummy",1000,3, 209,1000,0.11, 1.05);
  SetStyleHistoTH2ForGraphs(histResoDummy, "#it{E} (GeV)","< #it{E}_{rec} / #it{E}_{true} >", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.91);
  histResoDummy->GetXaxis()->SetNoExponent();
  histResoDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histResoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      histResoDummy->Draw();
      if(h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo]){
        DrawGammaSetMarker(h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo], 20, 3, colorClus[1][ialgo], colorClus[1][ialgo]);

        h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FHCALSTANDALONE_pion_meanE_E[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
      }
      if(h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo]){
        DrawGammaSetMarker(h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo], 20, 3, colorClus[0][ialgo], colorClus[0][ialgo]);

        h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FHCALSTANDALONE_pionComb_meanE_E[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
      }
      if(h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo]){
        DrawGammaSetMarker(h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo], 20, 3, colorClus[2][ialgo], colorClus[2][ialgo]);

        h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FHCALSTANDALONE_pionCombCalib_meanE_E[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
      }
      // cReso2x->Print(Form("%s/FHCALSTANDALONE_pion_MeanE_%s%s_Etrue.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
  }
}


  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=0)continue;
      cout << "PION!!!!" << endl;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      cout << str_clusterizer[ialgo].Data() << endl;
      histResoDummy->Draw();
      if(h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo]){
        DrawGammaSetMarker(h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo], 20, 3, colorClus[1][ialgo], colorClus[1][ialgo]);

        h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FHCALSTANDALONE_pion_meanE_Erec[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
      }
      if(h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo]){
        DrawGammaSetMarker(h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo], 20, 3, colorClus[2][ialgo], colorClus[2][ialgo]);

        h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FHCALSTANDALONE_pionComb_meanE_Erec[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
        for(int iparam=0;iparam<FunctionNL_PowExp150MeV->GetNpar();iparam++){
        cout << FunctionNL_PowExp150MeV->GetParameter(iparam) << ", ";
        }
        cout << endl;
      }
      if(h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo]){
        DrawGammaSetMarker(h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo], 20, 3, colorClus[2][ialgo], colorClus[2][ialgo]);

        h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo]->Draw("same,p");

        // Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        // TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        // FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        // h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        TF1* FunctionNLconst   = new TF1("FunctionNLconst",   "( [0] )",   1, 120);
        h_FHCALSTANDALONE_pionCombCalib_meanE_Erec[icalo][ialgo]->Fit(FunctionNLconst,"QNRME+","",1.0, 100);
        // FunctionNL_PowExp150MeV->Draw("same");
        FunctionNLconst->Draw("same");
      }
      cReso2x->Print(Form("%s/FHCALSTANDALONE_pion_MeanE_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
  }
}



  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=1)continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      histResoDummy->Draw();
      if(h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo]){
        DrawGammaSetMarker(h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo], 20, 3, colorClus[0][ialgo], colorClus[0][ialgo]);

        h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FEMCSTANDALONE_gamma_meanE_E[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
      }
      if(h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo]){
        DrawGammaSetMarker(h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo], 20, 3, colorClus[2][ialgo], colorClus[2][ialgo]);

        h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo]->Draw("same,p");

        // Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        // TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        // FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        // h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        TF1* FunctionNLconst   = new TF1("FunctionNLconst",   "( [0] )",   1, 120);
        h_FEMCSTANDALONE_gammaCalib_meanE_E[icalo][ialgo]->Fit(FunctionNLconst,"QNRME+","",1.0, 100);
        // FunctionNL_PowExp150MeV->Draw("same");
        FunctionNLconst->Draw("same");
      }
      // cReso2x->Print(Form("%s/FEMCSTANDALONE_gamma_MeanE_%s%s_Etrue.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
  }
}


  for(int icalo=0;icalo<_active_calo;icalo++){
    if(icalo!=1)continue;
      cout << "GAMMA!!!!" << endl;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(ialgo==kC3 || ialgo==k3x3) continue;
      cout << str_clusterizer[ialgo].Data() << endl;
      histResoDummy->Draw();
      if(h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo]){
        DrawGammaSetMarker(h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo], 20, 3, colorClus[0][ialgo], colorClus[0][ialgo]);

        h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FEMCSTANDALONE_gamma_meanE_Erec[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
        // for(int iparam=0;iparam<FunctionNL_PowExp150MeV->GetNpar();iparam++){
        // cout << FunctionNL_PowExp150MeV->GetParameter(iparam) << ", ";
        // }
        cout << endl;
      }
      if(h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo]){
        DrawGammaSetMarker(h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo], 20, 3, colorClus[2][ialgo], colorClus[2][ialgo]);

        h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo]->Draw("same,p");

        Double_t paramsTB150MeV[5]= {2.49398e+00,4.57472e-01,2.04677e+00,-6.90886e+02,9.11171e+02}; //noshift
        TF1* FunctionNL_PowExp150MeV   = new TF1("FunctionNL_PowExp150MeV",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
        FunctionNL_PowExp150MeV->SetParameters(paramsTB150MeV);
        h_FEMCSTANDALONE_gammaCalib_meanE_Erec[icalo][ialgo]->Fit(FunctionNL_PowExp150MeV,"QNRME+","",1.0, 100);
        FunctionNL_PowExp150MeV->Draw("same");
        for(int iparam=0;iparam<FunctionNL_PowExp150MeV->GetNpar();iparam++){
        cout << FunctionNL_PowExp150MeV->GetParameter(iparam) << ", ";
        }
        cout << endl;
      }
      cReso2x->Print(Form("%s/FEMCSTANDALONE_gamma_MeanE_%s%s_Erec.%s", outputDir.Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(), suffix.Data()));
  }
}

}

