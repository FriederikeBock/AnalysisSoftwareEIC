
// ANCHOR debug output verbosity
Int_t verbosityRH = 2;

// ANCHOR create histograms globally
TH2F*  h_RH_Reso_E[_active_calo][_active_algo];
TH2F*  h_RH_Reso_gamma_E[_active_calo][_active_algo];
TH2F*  h_RH_Reso_pion_E[_active_calo][_active_algo];
TH2F*  h_RH_Reso_gamma_E_Eta[_active_calo][_active_algo][nEta+1];
TH2F*  h_RH_Reso_pion_E_Eta[_active_calo][_active_algo][nEta+1];
TH2F*  h_RH_Reso_Erec[_active_calo][_active_algo];
TH2F*  h_RH_Reso_gamma_Erec[_active_calo][_active_algo];
TH2F*  h_RH_Reso_pion_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_gamma_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_pion_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_gamma_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_pion_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_gamma_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_pion_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_gamma_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_pion_E[_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_gamma_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_pion_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_gamma_Erec[_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_pion_Erec[_active_calo][_active_algo];

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

int nclusters_comb_RE = 0;
float* clusters_comb_RE_E            = new float[_maxNclusters];
float* clusters_comb_RE_Eta         = new float[_maxNclusters];
float* clusters_comb_RE_Phi         = new float[_maxNclusters];
float* clusters_comb_RE_M02         = new float[_maxNclusters];
float* clusters_comb_RE_M20         = new float[_maxNclusters];
bool* clusters_comb_RE_isMatched         = new bool[_maxNclusters];
int* clusters_comb_RE_NTower         = new int[_maxNclusters];
int* clusters_comb_RE_trueID       = new int[_maxNclusters];
int* clusters_comb_RE_NtrueID       = new int[_maxNclusters];
float* clusters_comb_RE_X         = new float[_maxNclusters];
float* clusters_comb_RE_Y         = new float[_maxNclusters];
float* clusters_comb_RE_Z         = new float[_maxNclusters];

// ANCHOR main function to be called in event loop
void resolutionhistos(){

  // f_FEMC_calib->SetParameters(paramsFEMCcalib);

  // f_FHCAL_calib->SetParameters(paramsHCALcalib);

  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
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
        int nbinsx = 400;
        if(icalo==kFHCAL) nbinsx = 200;
      if(!h_RH_Reso_E[icalo][ialgo])h_RH_Reso_E[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_Reso_gamma_E[icalo][ialgo])h_RH_Reso_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_Reso_pion_E[icalo][ialgo])h_RH_Reso_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      for (Int_t eT = 0; eT < nEta+1; eT++){
        if(!h_RH_Reso_gamma_E_Eta[icalo][ialgo][eT])h_RH_Reso_gamma_E_Eta[icalo][ialgo][eT] 	= new TH2F(Form("h_RH_Reso_gamma_E_Eta_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),eT), "", 500,0,250,100, 0, 2);
        if(!h_RH_Reso_pion_E_Eta[icalo][ialgo][eT])h_RH_Reso_pion_E_Eta[icalo][ialgo][eT] 	= new TH2F(Form("h_RH_Reso_pion_E_Eta_%s_%s_%d",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),eT), "", 500,0,250,100, 0, 2);
      }
      if(!h_RH_Reso_Erec[icalo][ialgo])h_RH_Reso_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_Reso_gamma_Erec[icalo][ialgo])h_RH_Reso_gamma_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_gamma_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_Reso_pion_Erec[icalo][ialgo])h_RH_Reso_pion_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_Reso_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);

      if(!h_RH_ResoCalib_E[icalo][ialgo])h_RH_ResoCalib_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_ResoCalib_gamma_E[icalo][ialgo])h_RH_ResoCalib_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_ResoCalib_pion_E[icalo][ialgo])h_RH_ResoCalib_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);

      if(!h_RH_ResoCalib_Erec[icalo][ialgo])h_RH_ResoCalib_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_ResoCalib_gamma_Erec[icalo][ialgo])h_RH_ResoCalib_gamma_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_gamma_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
      if(!h_RH_ResoCalib_pion_Erec[icalo][ialgo])h_RH_ResoCalib_pion_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCalib_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);


      for(Int_t iclus=0; iclus<nclusters_RH; iclus++){
        // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
        // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

        // cluster should have at least 2 towers
        if(clusters_RH_NTower[iclus]<2) continue;

        float clusterE_calo = clusters_RH_E[iclus];
        clusterE_calo/=getCalibrationValue(clusters_RH_E[iclus], icalo, ialgo);
        // loop over MC particles
        for(Int_t imc=0; imc<_nMCPart; imc++){

          // find true MC particle for given cluster
          if((clusters_RH_trueID[iclus])==_mcpart_ID[imc]){
            // if(verbosityRH>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;
            h_RH_Reso_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==22)h_RH_Reso_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==211)h_RH_Reso_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);

            Int_t et = 0;
            while ( ( clusters_RH_Eta[iclus] > partEta[et+1] ) && ( et < nEta )) et++;
            if(abs(_mcpart_PDG[imc])==22)h_RH_Reso_gamma_E_Eta[icalo][ialgo][et]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==211)h_RH_Reso_pion_E_Eta[icalo][ialgo][et]->Fill(_mcpart_E[imc],clusters_RH_E[iclus]/_mcpart_E[imc]);

            h_RH_Reso_Erec[icalo][ialgo]->Fill(clusters_RH_E[iclus],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==22)h_RH_Reso_gamma_Erec[icalo][ialgo]->Fill(clusters_RH_E[iclus],clusters_RH_E[iclus]/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==211)h_RH_Reso_pion_Erec[icalo][ialgo]->Fill(clusters_RH_E[iclus],clusters_RH_E[iclus]/_mcpart_E[imc]);


            h_RH_ResoCalib_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterE_calo/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==22)h_RH_ResoCalib_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterE_calo/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==211)h_RH_ResoCalib_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterE_calo/_mcpart_E[imc]);


            h_RH_ResoCalib_Erec[icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==22)h_RH_ResoCalib_gamma_Erec[icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[imc]);
            if(abs(_mcpart_PDG[imc])==211)h_RH_ResoCalib_pion_Erec[icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[imc]);
          }
        }
      }
    }
  }
  int icalo =0;
  for(int ialgo=0;ialgo<_active_algo;ialgo++){
    loadClusterizerInput(
      ialgo,
      icalo,
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
    loadClusterizerInput(
      ialgo,
      1,
      nclusters_comb_RE,
      clusters_comb_RE_E,
      clusters_comb_RE_Eta,
      clusters_comb_RE_Phi,
      clusters_comb_RE_M02,
      clusters_comb_RE_M20,
      clusters_comb_RE_isMatched,
      clusters_comb_RE_NTower,
      clusters_comb_RE_trueID,
      clusters_comb_RE_NtrueID,
      clusters_comb_RE_X,
      clusters_comb_RE_Y,
      clusters_comb_RE_Z);
    if(!h_RH_ResoComb_E[icalo][ialgo])h_RH_ResoComb_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoComb_gamma_E[icalo][ialgo])h_RH_ResoComb_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoComb_pion_E[icalo][ialgo])h_RH_ResoComb_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoCombCalib_E[icalo][ialgo])h_RH_ResoCombCalib_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoCombCalib_gamma_E[icalo][ialgo])h_RH_ResoCombCalib_gamma_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_gamma_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoCombCalib_pion_E[icalo][ialgo])h_RH_ResoCombCalib_pion_E[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_pion_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoComb_Erec[icalo][ialgo])h_RH_ResoComb_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoComb_gamma_Erec[icalo][ialgo])h_RH_ResoComb_gamma_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_gamma_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoComb_pion_Erec[icalo][ialgo])h_RH_ResoComb_pion_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoComb_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoCombCalib_Erec[icalo][ialgo])h_RH_ResoCombCalib_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoCombCalib_gamma_Erec[icalo][ialgo])h_RH_ResoCombCalib_gamma_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_gamma_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);
    if(!h_RH_ResoCombCalib_pion_Erec[icalo][ialgo])h_RH_ResoCombCalib_pion_Erec[icalo][ialgo] 	= new TH2F(Form("h_RH_ResoCombCalib_pion_Erec_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 500,0,250,400, 0, 2);


    for(Int_t iclus=0; iclus<nclusters_comb_RH; iclus++){
      // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
      // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

      // cluster should have at least 2 towers
      if(clusters_comb_RH_NTower[iclus]<2) continue;
      float clusterenergy = clusters_comb_RH_E[iclus];
      float clusterenergyCalib = clusters_comb_RH_E[iclus];
      clusterenergyCalib/=getCalibrationValue(clusters_comb_RH_E[iclus], kFHCAL, ialgo);
      // float clusterenergyCalib = clusters_comb_RH_E[iclus]/calibrationFHCAL(clusters_comb_RH_E[iclus]);;
      // loop over MC particles
      for(Int_t imc=0; imc<_nMCPart; imc++){

        // find true MC particle for given cluster
        if((clusters_comb_RH_trueID[iclus])==_mcpart_ID[imc]){

          for(Int_t iclus2=0; iclus2<nclusters_comb_RE; iclus2++){
            if(clusters_comb_RE_NTower[iclus2]<2) continue;
            if(clusters_comb_RH_trueID[iclus]==clusters_comb_RE_trueID[iclus2]){
              clusterenergy += (clusters_comb_RE_E[iclus2]);
              clusterenergyCalib += (clusters_comb_RE_E[iclus2]/=getCalibrationValue(clusters_comb_RE_E[iclus2], kFEMC, ialgo));
            }
          }
          // if(verbosityRH>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;

          h_RH_ResoComb_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergy/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==22)h_RH_ResoComb_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergy/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==211)h_RH_ResoComb_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergy/_mcpart_E[imc]);
          h_RH_ResoCombCalib_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergyCalib/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==22)h_RH_ResoCombCalib_gamma_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergyCalib/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==211)h_RH_ResoCombCalib_pion_E[icalo][ialgo]->Fill(_mcpart_E[imc],clusterenergyCalib/_mcpart_E[imc]);

          h_RH_ResoComb_Erec[icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==22)h_RH_ResoComb_gamma_Erec[icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==211)h_RH_ResoComb_pion_Erec[icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[imc]);
          h_RH_ResoCombCalib_Erec[icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==22)h_RH_ResoCombCalib_gamma_Erec[icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[imc]);
          if(abs(_mcpart_PDG[imc])==211)h_RH_ResoCombCalib_pion_Erec[icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[imc]);
        }
      }
    }
  }
}

// ANCHOR save function after event loop
void resolutionhistosSave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_RH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_RH_Reso_E[icalo][ialgo])h_RH_Reso_E[icalo][ialgo]->Write();
      if(h_RH_Reso_gamma_E[icalo][ialgo])h_RH_Reso_gamma_E[icalo][ialgo]->Write();
      if(h_RH_Reso_pion_E[icalo][ialgo])h_RH_Reso_pion_E[icalo][ialgo]->Write();

      for (Int_t eT = 0; eT < nEta+1; eT++){
        if(h_RH_Reso_gamma_E_Eta[icalo][ialgo][eT])h_RH_Reso_gamma_E_Eta[icalo][ialgo][eT]->Write();
        if(h_RH_Reso_pion_E_Eta[icalo][ialgo][eT])h_RH_Reso_pion_E_Eta[icalo][ialgo][eT]->Write();
      }

      if(h_RH_Reso_Erec[icalo][ialgo])h_RH_Reso_Erec[icalo][ialgo]->Write();
      if(h_RH_Reso_gamma_Erec[icalo][ialgo])h_RH_Reso_gamma_Erec[icalo][ialgo]->Write();
      if(h_RH_Reso_pion_Erec[icalo][ialgo])h_RH_Reso_pion_Erec[icalo][ialgo]->Write();

      if(h_RH_ResoCalib_E[icalo][ialgo])h_RH_ResoCalib_E[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_gamma_E[icalo][ialgo])h_RH_ResoCalib_gamma_E[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_pion_E[icalo][ialgo])h_RH_ResoCalib_pion_E[icalo][ialgo]->Write();

      if(h_RH_ResoCalib_Erec[icalo][ialgo])h_RH_ResoCalib_Erec[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_gamma_Erec[icalo][ialgo])h_RH_ResoCalib_gamma_Erec[icalo][ialgo]->Write();
      if(h_RH_ResoCalib_pion_Erec[icalo][ialgo])h_RH_ResoCalib_pion_Erec[icalo][ialgo]->Write();
      if(icalo==0){
        if(h_RH_ResoComb_E[icalo][ialgo])h_RH_ResoComb_E[icalo][ialgo]->Write();
        if(h_RH_ResoComb_gamma_E[icalo][ialgo])h_RH_ResoComb_gamma_E[icalo][ialgo]->Write();
        if(h_RH_ResoComb_pion_E[icalo][ialgo])h_RH_ResoComb_pion_E[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_E[icalo][ialgo])h_RH_ResoCombCalib_E[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_gamma_E[icalo][ialgo])h_RH_ResoCombCalib_gamma_E[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_pion_E[icalo][ialgo])h_RH_ResoCombCalib_pion_E[icalo][ialgo]->Write();
        if(h_RH_ResoComb_Erec[icalo][ialgo])h_RH_ResoComb_Erec[icalo][ialgo]->Write();
        if(h_RH_ResoComb_gamma_Erec[icalo][ialgo])h_RH_ResoComb_gamma_Erec[icalo][ialgo]->Write();
        if(h_RH_ResoComb_pion_Erec[icalo][ialgo])h_RH_ResoComb_pion_Erec[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_Erec[icalo][ialgo])h_RH_ResoCombCalib_Erec[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_gamma_Erec[icalo][ialgo])h_RH_ResoCombCalib_gamma_Erec[icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_pion_Erec[icalo][ialgo])h_RH_ResoCombCalib_pion_Erec[icalo][ialgo]->Write();
      }
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
