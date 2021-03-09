
// ANCHOR debug output verbosity
Int_t verbosityRH = 2;

// ANCHOR create histograms globally
TH2F*  h_RH_Reso_E[5][4];
TH2F*  h_RH_Reso_gamma_E[5][4];
TH2F*  h_RH_Reso_electron_E[5][4];
TH2F*  h_RH_Reso_pion_E[5][4];
TH2F*  h_RH_ResoCalib_E[5][4];
TH2F*  h_RH_ResoCalib_gamma_E[5][4];
TH2F*  h_RH_ResoCalib_electron_E[5][4];
TH2F*  h_RH_ResoCalib_pion_E[5][4];
TH2F*  h_RH_ResoComb_E[5][4];
TH2F*  h_RH_ResoComb_gamma_E[5][4];
TH2F*  h_RH_ResoComb_electron_E[5][4];
TH2F*  h_RH_ResoComb_pion_E[5][4];
TH2F*  h_RH_ResoCombCalib_E[5][4];
TH2F*  h_RH_ResoCombCalib_gamma_E[5][4];
TH2F*  h_RH_ResoCombCalib_electron_E[5][4];
TH2F*  h_RH_ResoCombCalib_pion_E[5][4];

// Double_t paramsFEMCcalib[5]= {3.15612e+00,1.57658e-01,1.66441e+00,-2.50979e+02,3.81511e+02}; //noshift
// TF1* f_FEMC_calib   = new TF1("f_FEMC_calib",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);
// Double_t paramsHCALcalib[5]= {1.75773e+00, 3.87923e-01, 1.20862e+00, -1.25194e+02, 7.54438e+01}; //noshift
// TF1* f_FHCAL_calib   = new TF1("f_FHCAL_calib",   "( [0] + [1] * TMath::Log(x) ) / ( 1 + ( [2] * TMath::Exp( ( x - [3] ) / [4] ) ) )",   1, 120);

// ANCHOR main function to be called in event loop
void resolutionhistos(){

  // f_FEMC_calib->SetParameters(paramsFEMCcalib);

  // f_FHCAL_calib->SetParameters(paramsHCALcalib);


  int nclusters_RH = 0;
  float* clusters_RH_E            = new float[_maxNclusters];
  float* clusters_RH_Eta         = new float[_maxNclusters];
  float* clusters_RH_Phi         = new float[_maxNclusters];
  float* clusters_RH_M02         = new float[_maxNclusters];
  float* clusters_RH_M20         = new float[_maxNclusters];
  bool* clusters_RH_isMatched         = new bool[_maxNclusters];
  int* clusters_RH_NTower         = new int[_maxNclusters];
  int* clusters_RH_trueID       = new int[_maxNclusters];
  int* clusters_RH_NtrueID       = new int[_maxNclusters];
  float* clusters_RH_X         = new float[_maxNclusters];
  float* clusters_RH_Y         = new float[_maxNclusters];
  float* clusters_RH_Z         = new float[_maxNclusters];
  for(int icalo=0;icalo<5;icalo++){
    for(int ialgo=0;ialgo<4;ialgo++){
      loadClusterizerInput(
        ialgo,
        icalo,
        nclusters_RH,
        clusters_RH_E,
        clusters_RH_Eta,
        clusters_RH_Phi,
        clusters_RH_M02,
        clusters_RH_M20,
        clusters_RH_isMatched,
        clusters_RH_NTower,
        clusters_RH_trueID,
        clusters_RH_NtrueID,
        clusters_RH_X,
        clusters_RH_Y,
        clusters_RH_Z);
      if(!h_RH_Reso_E[icalo][ialgo])h_RH_Reso_E[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
      if(!h_RH_Reso_gamma_E[icalo][ialgo])h_RH_Reso_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
      if(!h_RH_Reso_electron_E[icalo][ialgo])h_RH_Reso_electron_E[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_electron_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
      if(!h_RH_Reso_pion_E[icalo][ialgo])h_RH_Reso_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
      if(!h_RH_ResoCalib_E[icalo][ialgo])h_RH_ResoCalib_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
      if(!h_RH_ResoCalib_gamma_E[icalo][ialgo])h_RH_ResoCalib_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
      if(!h_RH_ResoCalib_electron_E[icalo][ialgo])h_RH_ResoCalib_electron_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_electron_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
      if(!h_RH_ResoCalib_pion_E[icalo][ialgo])h_RH_ResoCalib_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);


      for(Int_t iclus=0; iclus<nclusters_RH; iclus++){
        // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
        // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

        // cluster should have at least 2 towers
        if(clusters_RH_NTower[iclus]<2) continue;

        // loop over MC particles
        for(Int_t imc=0; imc<_nMCPart; imc++){

          // find true MC particle for given cluster
          if((clusters_RH_trueID[iclus])==_mcpart_ID[imc]){
            // if(verbosityRH>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;
            h_RH_Reso_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==22)h_RH_Reso_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==11)h_RH_Reso_electron_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==211)h_RH_Reso_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);


            float clusterE_calo = clusters_RH_E[iclus];
            // if(icalo==kFHCAL)
            //   clusterE_calo/=calibrationFHCAL(clusters_RH_E[iclus]);
            // else if(icalo==kFEMC)
            //   clusterE_calo/=calibrationFEMC(clusters_RH_E[iclus]);

            // if(icalo==kFHCAL)clusterE_calo/=f_FHCAL_calib->Eval(clusters_RH_E[iclus]);
            // if(icalo==kFEMC)clusterE_calo/=f_FEMC_calib->Eval(clusters_RH_E[iclus]);
            h_RH_ResoCalib_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterE_calo/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==22)h_RH_ResoCalib_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterE_calo/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==11)h_RH_ResoCalib_electron_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterE_calo/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==211)h_RH_ResoCalib_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterE_calo/_mcpart_E[imc]);
          }
        }
      }
    }
  }
  int icalo = 0;
  int nclusters_comb_RH = 0;
  float* clusters_comb_RH_E            = new float[_maxNclusters];
  float* clusters_comb_RH_Eta         = new float[_maxNclusters];
  float* clusters_comb_RH_Phi         = new float[_maxNclusters];
  float* clusters_comb_RH_M02         = new float[_maxNclusters];
  float* clusters_comb_RH_M20         = new float[_maxNclusters];
  bool* clusters_comb_RH_isMatched         = new bool[_maxNclusters];
  int* clusters_comb_RH_NTower         = new int[_maxNclusters];
  int* clusters_comb_RH_trueID       = new int[_maxNclusters];
  int* clusters_comb_RH_NtrueID       = new int[_maxNclusters];
  float* clusters_comb_RH_X         = new float[_maxNclusters];
  float* clusters_comb_RH_Y         = new float[_maxNclusters];
  float* clusters_comb_RH_Z         = new float[_maxNclusters];
  for(int ialgo=0;ialgo<4;ialgo++){
    loadClusterizerInput(
      ialgo,
      icalo,
      nclusters_RH,
      clusters_RH_E,
      clusters_RH_Eta,
      clusters_RH_Phi,
      clusters_RH_M02,
      clusters_RH_M20,
      clusters_RH_isMatched,
      clusters_RH_NTower,
      clusters_RH_trueID,
      clusters_RH_NtrueID,
      clusters_RH_X,
      clusters_RH_Y,
      clusters_RH_Z);
    loadClusterizerInput(
      ialgo,
      1,
      nclusters_comb_RH,
      clusters_comb_RH_E,
      clusters_comb_RH_Eta,
      clusters_comb_RH_Phi,
      clusters_comb_RH_M02,
      clusters_comb_RH_M20,
      clusters_comb_RH_isMatched,
      clusters_comb_RH_NTower,
      clusters_comb_RH_trueID,
      clusters_comb_RH_NtrueID,
      clusters_comb_RH_X,
      clusters_comb_RH_Y,
      clusters_comb_RH_Z);
    if(!h_RH_ResoComb_E[icalo][ialgo])h_RH_ResoComb_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
    if(!h_RH_ResoComb_gamma_E[icalo][ialgo])h_RH_ResoComb_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
    if(!h_RH_ResoComb_electron_E[icalo][ialgo])h_RH_ResoComb_electron_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_electron_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
    if(!h_RH_ResoComb_pion_E[icalo][ialgo])h_RH_ResoComb_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
    if(!h_RH_ResoCombCalib_E[icalo][ialgo])h_RH_ResoCombCalib_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
    if(!h_RH_ResoCombCalib_gamma_E[icalo][ialgo])h_RH_ResoCombCalib_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
    if(!h_RH_ResoCombCalib_electron_E[icalo][ialgo])h_RH_ResoCombCalib_electron_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_electron_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);
    if(!h_RH_ResoCombCalib_pion_E[icalo][ialgo])h_RH_ResoCombCalib_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 240,0,120, 100, 0, 2);


    for(Int_t iclus=0; iclus<nclusters_RH; iclus++){
      // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
      // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

      // cluster should have at least 2 towers
      if(clusters_RH_NTower[iclus]<2) continue;
      float clusterenergy = clusters_RH_E[iclus];
      float clusterenergyCalib = clusters_RH_E[iclus];
      // float clusterenergyCalib = clusters_RH_E[iclus]/calibrationFHCAL(clusters_RH_E[iclus]);;
      // loop over MC particles
      for(Int_t imc=0; imc<_nMCPart; imc++){

        // find true MC particle for given cluster
        if((clusters_RH_trueID[iclus])==_mcpart_ID[imc]){

          for(Int_t iclus2=0; iclus2<nclusters_comb_RH; iclus2++){
            if(clusters_comb_RH_NTower[iclus]<2) continue;
            if(clusters_RH_trueID[iclus]==clusters_comb_RH_trueID[iclus]){
              clusterenergy += (clusters_comb_RH_E[iclus]);
              clusterenergyCalib += (clusters_comb_RH_E[iclus]);
              // clusterenergy += (clusters_comb_RH_E[iclus]/calibrationFEMC(clusters_comb_RH_E[iclus]));
              // clusterenergyCalib += (clusters_comb_RH_E[iclus]/calibrationFEMC(clusters_comb_RH_E[iclus]));
            }
          }
          // if(verbosityRH>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;

          h_RH_ResoComb_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergy/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==22)h_RH_ResoComb_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergy/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==11)h_RH_ResoComb_electron_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergy/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==211)h_RH_ResoComb_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergy/_mcpart_E[imc]);
          h_RH_ResoCombCalib_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergyCalib/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==22)h_RH_ResoCombCalib_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergyCalib/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==11)h_RH_ResoCombCalib_electron_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergyCalib/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==211)h_RH_ResoCombCalib_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergyCalib/_mcpart_E[imc]);
        }
      }
    }
  }
}

// ANCHOR save function after event loop
void resolutionhistosSave(){
  // make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/resolutionhistos");

  // gSystem->Exec("mkdir -p treeProcessing/resolutionhistos");
  // define output file
  TFile* fileOutput = new TFile(Form("%s/resolutionhistos/output_RH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<5;icalo++){
    for(int ialgo=0;ialgo<4;ialgo++){
      if(h_RH_Reso_E[icalo][ialgo])h_RH_Reso_E[icalo][ialgo]->Write();
      if(h_RH_Reso_gamma_E[icalo][ialgo])h_RH_Reso_gamma_E[icalo][ialgo]->Write();
      if(h_RH_Reso_electron_E[icalo][ialgo])h_RH_Reso_electron_E[icalo][ialgo]->Write();
      if(h_RH_Reso_pion_E[icalo][ialgo])h_RH_Reso_pion_E[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_E[icalo][ialgo])h_RH_ResoCalib_E[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_gamma_E[icalo][ialgo])h_RH_ResoCalib_gamma_E[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_electron_E[icalo][ialgo])h_RH_ResoCalib_electron_E[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_pion_E[icalo][ialgo])h_RH_ResoCalib_pion_E[icalo][ialgo]->Write();
      if(icalo==0){
        if(h_RH_ResoComb_E[icalo][ialgo])h_RH_ResoComb_E[icalo][ialgo]->Write();
        if(h_RH_ResoComb_gamma_E[icalo][ialgo])h_RH_ResoComb_gamma_E[icalo][ialgo]->Write();
        if(h_RH_ResoComb_electron_E[icalo][ialgo])h_RH_ResoComb_electron_E[icalo][ialgo]->Write();
        if(h_RH_ResoComb_pion_E[icalo][ialgo])h_RH_ResoComb_pion_E[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_E[icalo][ialgo])h_RH_ResoCombCalib_E[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_gamma_E[icalo][ialgo])h_RH_ResoCombCalib_gamma_E[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_electron_E[icalo][ialgo])h_RH_ResoCombCalib_electron_E[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_pion_E[icalo][ialgo])h_RH_ResoCombCalib_pion_E[icalo][ialgo]->Write();
      }
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}