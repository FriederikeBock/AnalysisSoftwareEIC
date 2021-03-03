
// ANCHOR debug output verbosity
int int_CS_verbosity = 2;

int int_CS_mcparticles_PDG[6] = {11 /*e*/,22 /*gamma*/,211  /*pi*/,  2212/*p*/,  2112/*n*/,  321/*K*/};

// ANCHOR create histograms globally
TH1F*  h_CS_cluster_fhcal_V1_E 	= new TH1F("h_CS_cluster_fhcal_V1_E", "", 120, 0,60);
TH1F*  h_CS_cluster_fhcal_NxN_E 	= new TH1F("h_CS_cluster_fhcal_NxN_E", "", 120, 0,60);
TH1F*  h_CS_cluster_fhcal_XN_E 	= new TH1F("h_CS_cluster_fhcal_XN_E", "", 120, 0,60);
TH1F*  h_CS_cluster_fhcal_V3_E 	= new TH1F("h_CS_cluster_fhcal_V3_E", "", 120, 0,60);
TH2F*  h_CS_cluster_M02_fhcal_V1_E 	= new TH2F("h_CS_cluster_M02_fhcal_V1_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_M02_fhcal_NxN_E 	= new TH2F("h_CS_cluster_M02_fhcal_NxN_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_M02_fhcal_XN_E 	= new TH2F("h_CS_cluster_M02_fhcal_XN_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_M02_fhcal_V3_E 	= new TH2F("h_CS_cluster_M02_fhcal_V3_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_M02_fhcal_V1_E_truth[6] = {NULL};
TH2F*  h_CS_cluster_M02_fhcal_V3_E_truth[6] = {NULL};
TH2F*  h_CS_cluster_M02_fhcal_NxN_E_truth[6] = {NULL};
TH2F*  h_CS_cluster_M02_fhcal_XN_E_truth[6] = {NULL};
TH2F*  h_CS_cluster_M20_fhcal_V1_E 	= new TH2F("h_CS_cluster_M20_fhcal_V1_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_M20_fhcal_NxN_E 	= new TH2F("h_CS_cluster_M20_fhcal_NxN_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_M20_fhcal_XN_E 	= new TH2F("h_CS_cluster_M20_fhcal_XN_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_M20_fhcal_V3_E 	= new TH2F("h_CS_cluster_M20_fhcal_V3_E", "", 100, 0, 2, 120, 0,60);
TH2F*  h_CS_cluster_NTower_fhcal_V1_E 	= new TH2F("h_CS_cluster_nTower_fhcal_V1_E", "", 120, 0,60, 30, -0.5, 29.5);
TH2F*  h_CS_cluster_nTower_fhcal_NxN_E 	= new TH2F("h_CS_cluster_nTower_fhcal_NxN_E", "", 120, 0,60, 30, -0.5, 29.5);
TH2F*  h_CS_cluster_nTower_fhcal_XN_E 	= new TH2F("h_CS_cluster_nTower_fhcal_XN_E", "", 120, 0,60, 30, -0.5, 29.5);
TH2F*  h_CS_cluster_nTower_fhcal_V3_E 	= new TH2F("h_CS_cluster_nTower_fhcal_V3_E", "", 120, 0,60, 30, -0.5, 29.5);

// ANCHOR main function to be called in event loop
void clusterstudies(){

  for(Int_t iclus=0; iclus<_nclusters_FHCAL; iclus++){
    // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _cluster_FHCAL_E[iclus] << " GeV" << endl;
    // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster MC ID " << _cluster_FHCAL_trueID[iclus] << endl;

    // cluster should have at least 2 towers
    // if(_cluster_FHCAL_NTower[iclus]<2) continue;
    h_CS_cluster_NTower_fhcal_V1_E->Fill(_cluster_FHCAL_E[iclus],_cluster_FHCAL_NTower[iclus]);

    // // loop over MC particles
    // for(Int_t imc=0; imc<_nMCPart; imc++){

    //   // find true MC particle for given cluster
    //   if(_cluster_FHCAL_trueID[iclus]==_mcpart_ID[imc]){
    //     // if(int_CS_verbosity>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;

    //     // h_CS_cluster_M02_fhcal_V1_E->Fill(_mcpart_E[imc],_cluster_FHCAL_E[iclus]/_mcpart_E[imc]);
    //   }
    // }
  }
  if(_do_NxNclusterizer){
    for(Int_t iclus=0; iclus<_nclusters_NxN_FHCAL; iclus++){
      // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _cluster_FHCAL_E[iclus] << " GeV" << endl;
      // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster MC ID " << _cluster_FHCAL_trueID[iclus] << endl;

      h_CS_cluster_fhcal_NxN_E->Fill(_clusters_NxN_FHCAL_E[iclus]);
      h_CS_cluster_M02_fhcal_NxN_E->Fill(_clusters_NxN_FHCAL_M02[iclus],_clusters_NxN_FHCAL_E[iclus]);
      h_CS_cluster_M20_fhcal_NxN_E->Fill(_clusters_NxN_FHCAL_M20[iclus],_clusters_NxN_FHCAL_E[iclus]);

      h_CS_cluster_nTower_fhcal_NxN_E->Fill(_clusters_NxN_FHCAL_E[iclus],_clusters_NxN_FHCAL_NTower[iclus]);

      for (int iPDG = 0; iPDG < 6; iPDG++)
      {
        if(!h_CS_cluster_M02_fhcal_NxN_E_truth[iPDG]) h_CS_cluster_M02_fhcal_NxN_E_truth[iPDG] = new TH2F(Form("h_CS_cluster_M02_fhcal_NxN_E_truth_%d",int_CS_mcparticles_PDG[iPDG]), "", 100, 0, 2, 120, 0,60);
        if(abs(_mcpart_PDG[_clusters_NxN_FHCAL_trueID[iclus]-1]) == abs(int_CS_mcparticles_PDG[iPDG])){
          h_CS_cluster_M02_fhcal_NxN_E_truth[iPDG]->Fill(_clusters_NxN_FHCAL_M02[iclus],_clusters_NxN_FHCAL_E[iclus]);
        }
      }
    }
  }
  if(_do_XNclusterizer){
    for(Int_t iclus=0; iclus<_nclusters_XN_FHCAL; iclus++){
      // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _cluster_FHCAL_E[iclus] << " GeV" << endl;
      // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster MC ID " << _cluster_FHCAL_trueID[iclus] << endl;

      h_CS_cluster_fhcal_XN_E->Fill(_clusters_XN_FHCAL_E[iclus]);
      h_CS_cluster_M02_fhcal_XN_E->Fill(_clusters_XN_FHCAL_M02[iclus],_clusters_XN_FHCAL_E[iclus]);
      h_CS_cluster_M20_fhcal_XN_E->Fill(_clusters_XN_FHCAL_M20[iclus],_clusters_XN_FHCAL_E[iclus]);

      h_CS_cluster_nTower_fhcal_XN_E->Fill(_clusters_XN_FHCAL_E[iclus],_clusters_XN_FHCAL_NTower[iclus]);

      for (int iPDG = 0; iPDG < 6; iPDG++)
      {
        if(!h_CS_cluster_M02_fhcal_XN_E_truth[iPDG]) h_CS_cluster_M02_fhcal_XN_E_truth[iPDG] = new TH2F(Form("h_CS_cluster_M02_fhcal_XN_E_truth_%d",int_CS_mcparticles_PDG[iPDG]), "", 100, 0, 2, 120, 0,60);
        if(abs(_mcpart_PDG[_clusters_XN_FHCAL_trueID[iclus]-1]) == abs(int_CS_mcparticles_PDG[iPDG])){
          h_CS_cluster_M02_fhcal_XN_E_truth[iPDG]->Fill(_clusters_XN_FHCAL_M02[iclus],_clusters_XN_FHCAL_E[iclus]);
        }
      }
    }
  }
  if(_do_V3clusterizer){
    for(Int_t iclus=0; iclus<_nclusters_V3_FHCAL; iclus++){
      // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _cluster_FHCAL_E[iclus] << " GeV" << endl;
      // if(int_CS_verbosity>1) cout << "\tFHCAL: cluster MC ID " << _cluster_FHCAL_trueID[iclus] << endl;

      h_CS_cluster_fhcal_V3_E->Fill(_clusters_V3_FHCAL_E[iclus]);

      h_CS_cluster_M02_fhcal_V3_E->Fill(_clusters_V3_FHCAL_M02[iclus],_clusters_V3_FHCAL_E[iclus]);
      h_CS_cluster_M20_fhcal_V3_E->Fill(_clusters_V3_FHCAL_M20[iclus],_clusters_V3_FHCAL_E[iclus]);

      h_CS_cluster_nTower_fhcal_V3_E->Fill(_clusters_V3_FHCAL_E[iclus],_clusters_V3_FHCAL_NTower[iclus]);

      // cluster should have at least 2 towers
      // if(_clusters_V3_FHCAL_NTower[iclus]<2) continue;
      // if(_clusters_V3_FHCAL_M02[iclus]>1)cout << "clsTrueID: " << _clusters_V3_FHCAL_trueID[iclus] << "\tnIDMC: " << _nMCPart<< "\tPDG: "   << abs(_mcpart_PDG[_clusters_V3_FHCAL_trueID[iclus]]) << endl;
      // find true MC particle for given cluster
      for (int iPDG = 0; iPDG < 6; iPDG++)
      {
        if(!h_CS_cluster_M02_fhcal_V3_E_truth[iPDG]) h_CS_cluster_M02_fhcal_V3_E_truth[iPDG] = new TH2F(Form("h_CS_cluster_M02_fhcal_V3_E_truth_%d",int_CS_mcparticles_PDG[iPDG]), "", 100, 0, 2, 120, 0,60);
        if(abs(_mcpart_PDG[_clusters_V3_FHCAL_trueID[iclus]-1]) == abs(int_CS_mcparticles_PDG[iPDG])){
          h_CS_cluster_M02_fhcal_V3_E_truth[iPDG]->Fill(_clusters_V3_FHCAL_M02[iclus],_clusters_V3_FHCAL_E[iclus]);
        }
      }
    }
  }


}

// ANCHOR save function after event loop
void clusterstudiesSave(){
  // make output directory
    gSystem->Exec("mkdir -p "+outputDir + "/clusterstudies");

  // gSystem->Exec("mkdir -p treeProcessing/clusterstudies");
  // define output file
  TFile* fileOutput = new TFile(Form("%s/clusterstudies/output_RH.root",outputDir.Data()),"RECREATE");

  // write histograms
  if(h_CS_cluster_fhcal_NxN_E) h_CS_cluster_fhcal_NxN_E->Write();
  if(h_CS_cluster_fhcal_XN_E) h_CS_cluster_fhcal_XN_E->Write();
  if(h_CS_cluster_fhcal_V3_E) h_CS_cluster_fhcal_V3_E->Write();
  if(h_CS_cluster_NTower_fhcal_V1_E) h_CS_cluster_NTower_fhcal_V1_E->Write();
  if(h_CS_cluster_nTower_fhcal_NxN_E) h_CS_cluster_nTower_fhcal_NxN_E->Write();
  if(h_CS_cluster_nTower_fhcal_XN_E) h_CS_cluster_nTower_fhcal_XN_E->Write();
  if(h_CS_cluster_nTower_fhcal_V3_E) h_CS_cluster_nTower_fhcal_V3_E->Write();

  if(h_CS_cluster_M02_fhcal_V1_E) h_CS_cluster_M02_fhcal_V1_E->Write();
  if(h_CS_cluster_M02_fhcal_NxN_E) h_CS_cluster_M02_fhcal_NxN_E->Write();
  if(h_CS_cluster_M02_fhcal_XN_E) h_CS_cluster_M02_fhcal_XN_E->Write();
  if(h_CS_cluster_M02_fhcal_V3_E) h_CS_cluster_M02_fhcal_V3_E->Write();
  for (int iPDG = 0; iPDG < 6; iPDG++)
  {
      if(h_CS_cluster_M02_fhcal_V3_E_truth[iPDG]) h_CS_cluster_M02_fhcal_V3_E_truth[iPDG]->Write();;
      if(h_CS_cluster_M02_fhcal_NxN_E_truth[iPDG]) h_CS_cluster_M02_fhcal_NxN_E_truth[iPDG]->Write();;
      if(h_CS_cluster_M02_fhcal_XN_E_truth[iPDG]) h_CS_cluster_M02_fhcal_XN_E_truth[iPDG]->Write();;
  }
  if(h_CS_cluster_M20_fhcal_V1_E) h_CS_cluster_M20_fhcal_V1_E->Write();
  if(h_CS_cluster_M20_fhcal_NxN_E) h_CS_cluster_M20_fhcal_NxN_E->Write();
  if(h_CS_cluster_M20_fhcal_XN_E) h_CS_cluster_M20_fhcal_XN_E->Write();
  if(h_CS_cluster_M20_fhcal_V3_E) h_CS_cluster_M20_fhcal_V3_E->Write();

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}