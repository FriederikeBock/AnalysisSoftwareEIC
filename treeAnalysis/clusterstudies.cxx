
// ANCHOR debug output verbosity
int int_CS_verbosity = 2;

const int nParticlesPDG_CS = 6;
int int_CS_mcparticles_PDG[nParticlesPDG_CS] = {11 /*e*/,22 /*gamma*/,211  /*pi*/,  2212/*p*/,  2112/*n*/,  321/*K*/};
TString str_CS_mcparticles[nParticlesPDG_CS] = {"electron", "photon", "cpion", "proton", "neutron", "kaon"};

// ANCHOR create histograms globally
TH1F*  h_CS_clusters_E[_active_calo][_active_algo] = {{NULL}};
TH2F*  h_CS_clusters_M02_E[_active_calo][_active_algo] = {{NULL}};
TH2F*  h_CS_clusters_M02_part_E_truth[_active_calo][_active_algo][nParticlesPDG_CS] = {{{NULL}}};
TH2F*  h_CS_clusters_M20_E[_active_calo][_active_algo] = {{NULL}};
TH2F*  h_CS_clusters_M20_part_E_truth[_active_calo][_active_algo][nParticlesPDG_CS] = {{{NULL}}};
TH2F*  h_CS_clusters_NTower_E[_active_calo][_active_algo] = {{NULL}};

int nclusters_CLSTD = 0;
float* clusters_CLSTD_E         = new float[_maxNclusters];
float* clusters_CLSTD_Eta       = new float[_maxNclusters];
float* clusters_CLSTD_Phi       = new float[_maxNclusters];
float* clusters_CLSTD_M02       = new float[_maxNclusters];
float* clusters_CLSTD_M20       = new float[_maxNclusters];
bool* clusters_CLSTD_isMatched  = new bool[_maxNclusters];
int* clusters_CLSTD_NTower      = new int[_maxNclusters];
int* clusters_CLSTD_trueID      = new int[_maxNclusters];
int* clusters_CLSTD_NtrueID     = new int[_maxNclusters];
float* clusters_CLSTD_X         = new float[_maxNclusters];
float* clusters_CLSTD_Y         = new float[_maxNclusters];
float* clusters_CLSTD_Z         = new float[_maxNclusters];

// ANCHOR main function to be called in event loop
void clusterstudies(){
  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){

      // initialize histos if not yet done
      if(!h_CS_clusters_E[icalo][ialgo])h_CS_clusters_E[icalo][ialgo]= new TH1F(Form("h_CS_clusters_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP);
      if(!h_CS_clusters_M02_E[icalo][ialgo])h_CS_clusters_M02_E[icalo][ialgo]= new TH2F(Form("h_CS_clusters_M02_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "" , 100, 0, 2, nBinsP, binningP);
      if(!h_CS_clusters_M20_E[icalo][ialgo])h_CS_clusters_M20_E[icalo][ialgo]= new TH2F(Form("h_CS_clusters_M20_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "" , 100, 0, 2, nBinsP, binningP);
      if(!h_CS_clusters_NTower_E[icalo][ialgo])h_CS_clusters_NTower_E[icalo][ialgo]= new TH2F(Form("h_CS_clusters_NTower_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "" , nBinsP, binningP, 30, -0.5, 29.5);
      for (int iPDG = 0; iPDG < nParticlesPDG_CS; iPDG++){
        if(!h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG])h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]= new TH2F(Form("h_CS_clusters_M02_part_E_truth_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "" , 100, 0, 2, nBinsP, binningP);
        if(!h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG])h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]= new TH2F(Form("h_CS_clusters_M20_part_E_truth_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "" , 100, 0, 2, nBinsP, binningP);
      }

      loadClusterizerInput(
        ialgo,
        icalo,
        nclusters_CLSTD,
        clusters_CLSTD_E,
        clusters_CLSTD_Eta,
        clusters_CLSTD_Phi,
        clusters_CLSTD_M02,
        clusters_CLSTD_M20,
        clusters_CLSTD_isMatched,
        clusters_CLSTD_NTower,
        clusters_CLSTD_trueID,
        clusters_CLSTD_NtrueID,
        clusters_CLSTD_X,
        clusters_CLSTD_Y,
        clusters_CLSTD_Z
      );

      for(Int_t iclus=0; iclus<nclusters_CLSTD; iclus++){
        h_CS_clusters_NTower_E[icalo][ialgo]->Fill(clusters_CLSTD_E[iclus],clusters_CLSTD_NTower[iclus]);

        h_CS_clusters_E[icalo][ialgo]->Fill(clusters_CLSTD_E[iclus]);


        h_CS_clusters_M02_E[icalo][ialgo]->Fill(clusters_CLSTD_M02[iclus],clusters_CLSTD_E[iclus]);
        h_CS_clusters_M20_E[icalo][ialgo]->Fill(clusters_CLSTD_M20[iclus],clusters_CLSTD_E[iclus]);

        for (int iPDG = 0; iPDG < nParticlesPDG_CS; iPDG++)
        {
          if(abs(_mcpart_PDG[clusters_CLSTD_NtrueID[iclus]-1]) == abs(int_CS_mcparticles_PDG[iPDG])){
            h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]->Fill(clusters_CLSTD_M02[iclus],clusters_CLSTD_E[iclus]);
            h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]->Fill(clusters_CLSTD_M20[iclus],clusters_CLSTD_E[iclus]);
          }
        }
      }
    }
  }
}

// ANCHOR save function after event loop
void clusterstudiesSave(){
  // gSystem->Exec("mkdir -p treeProcessing/clusterstudies");
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_RH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<_active_calo;icalo++){
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_CS_clusters_NTower_E[icalo][ialgo]) h_CS_clusters_NTower_E[icalo][ialgo]->Write();
      if(h_CS_clusters_E[icalo][ialgo]) h_CS_clusters_E[icalo][ialgo]->Write();
      if(h_CS_clusters_M02_E[icalo][ialgo]) h_CS_clusters_M02_E[icalo][ialgo]->Write();
      if(h_CS_clusters_M20_E[icalo][ialgo]) h_CS_clusters_M20_E[icalo][ialgo]->Write();

      for (int iPDG = 0; iPDG < nParticlesPDG_CS; iPDG++){
        if(h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]) h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]->Write();
        if(h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]) h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]->Write();
      }
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
