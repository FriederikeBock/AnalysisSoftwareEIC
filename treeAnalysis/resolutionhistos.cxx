
// ANCHOR debug output verbosity
Int_t verbosityRH = 2;

// ANCHOR create histograms globally
TH2F*  h_RH_pion_fhcal_V1_E 	= new TH2F("h_RH_pion_fhcal_V1_E", "", 120, 0,60, 100, 0, 2);
TH2F*  h_RH_pion_fhcal_NxN_E 	= new TH2F("h_RH_pion_fhcal_NxN_E", "", 120, 0,60, 100, 0, 2);
TH2F*  h_RH_pion_fhcal_V3_E 	= new TH2F("h_RH_pion_fhcal_V3_E", "", 120, 0,60, 100, 0, 2);

// ANCHOR main function to be called in event loop
void resolutionhistos(){


  for(Int_t iclus=0; iclus<_nclusters_FHCAL; iclus++){
    // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
    // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

    // cluster should have at least 2 towers
    if(_clusters_FHCAL_NTower[iclus]<2) continue;

    // loop over MC particles
    for(Int_t imc=0; imc<_nMCPart; imc++){

      // find true MC particle for given cluster
      if((_clusters_FHCAL_trueID[iclus]-1)==_mcpart_ID[imc]){
        // if(verbosityRH>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;

        h_RH_pion_fhcal_V1_E->Fill(_mcpart_E[imc],_clusters_FHCAL_E[iclus]/_mcpart_E[imc]);
      }
    }
  }
  if(_do_NxNclusterizer){
    for(Int_t iclus=0; iclus<_nclusters_NxN_FHCAL; iclus++){
      // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
      // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

      // cluster should have at least 2 towers
      if(_clusters_NxN_FHCAL_NTower[iclus]<2) continue;

      // loop over MC particles
      for(Int_t imc=0; imc<_nMCPart; imc++){

        // find true MC particle for given cluster
        if((_clusters_NxN_FHCAL_trueID[iclus]-1)==_mcpart_ID[imc]){
          // if(verbosityRH>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;

          h_RH_pion_fhcal_NxN_E->Fill(_mcpart_E[imc],_clusters_NxN_FHCAL_E[iclus]/_mcpart_E[imc]);
        }
      }
    }
  }
  if(_do_V3clusterizer){
    for(Int_t iclus=0; iclus<_nclusters_V3_FHCAL; iclus++){
      // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
      // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

      // cluster should have at least 2 towers
      if(_clusters_V3_FHCAL_NTower[iclus]<2) continue;

      // loop over MC particles
      for(Int_t imc=0; imc<_nMCPart; imc++){

        // find true MC particle for given cluster
        if((_clusters_V3_FHCAL_trueID[iclus]-1)==_mcpart_ID[imc]){
          // if(verbosityRH>1) cout << "\tfound MC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;

          h_RH_pion_fhcal_V3_E->Fill(_mcpart_E[imc],_clusters_V3_FHCAL_E[iclus]/_mcpart_E[imc]);
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
  h_RH_pion_fhcal_V1_E->Write();
  h_RH_pion_fhcal_NxN_E->Write();
  h_RH_pion_fhcal_V3_E->Write();

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}