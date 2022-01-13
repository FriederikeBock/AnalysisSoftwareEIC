TFile* epidFile = nullptr;
TTree* epidTree = nullptr;

double p_track;
double eta_track;
double pT_track;
double theta_track;
double phi_track;
double E_clusters_all[maxcalo];	
double E_clusters_match[maxcalo];	
double Emax_cluster_M02;
double Emax_cluster_M20;
double E_MC;

void reset_variables();

void make_epid_tree(TString epidFilename) {
  
  epidFile = new TFile(epidFilename, "RECREATE");
  epidTree = new TTree("T", "Track and cluster tree for ePID");
  epidTree->Branch("p_track",		&p_track,		"p_track/D");
  epidTree->Branch("eta_track",		&eta_track,		"eta_track/D");
  epidTree->Branch("pT_track",		&pT_track,		"pT_track/D");
  epidTree->Branch("theta_track",		&theta_track,		"theta_track/D");
  epidTree->Branch("phi_track",		&phi_track,		"phi_track/D");
  epidTree->Branch("E_clusters_all",	E_clusters_all,		"E_clusters_all[11]/D");
  epidTree->Branch("E_clusters_match",	E_clusters_match,	"E_clusters_match[11]/D");
  epidTree->Branch("Emax_cluster_M02",	&Emax_cluster_M02,	"Emax_cluster_M02/D");	
  epidTree->Branch("Emax_cluster_M20",	&Emax_cluster_M20,	"Emax_cluster_M20/D");	
  epidTree->Branch("E_MC",		&E_MC,			"E_MC/D");
}

void electronpid( ) {

  Int_t fTrackSource = 0;
  Int_t fClusterAlgo = 0;

  reset_variables();

  for(Int_t itrk = 0; itrk < _nTracks; itrk++) {

    if(_track_source[itrk] != fTrackSource) continue;

    // Get track info
    TVector3 track_vec(_track_px[itrk], _track_py[itrk], _track_pz[itrk]);
    p_track = track_vec.Mag();
    eta_track = track_vec.PseudoRapidity();
    pT_track = track_vec.Pt();
    theta_track = track_vec.Theta();
    phi_track = track_vec.Phi();

    double maxClusterE = 0.;

    if (_track_matchECal[itrk].caloid == -1) continue;
    if (_track_matchECal[itrk].id == -1) continue;

    for (Int_t icalo = 0; icalo < maxcalo; icalo++) {		
      
      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[fClusterAlgo][icalo].size(); iclus++) {

        E_clusters_all[icalo]+=(_clusters_calo[fClusterAlgo][icalo].at(iclus)).cluster_E;

        if ((_clusters_calo[fClusterAlgo][icalo].at(iclus)).cluster_isMatched) {

          for(Int_t itrm = 0; itrm < (int)(_clusters_calo[fClusterAlgo][icalo].at(iclus)).cluster_matchedTracks.size(); itrm++) {
            int trkIDCl = ((_clusters_calo[fClusterAlgo][icalo].at(iclus)).cluster_matchedTracks).at(itrm).id;            
            if( trkIDCl == (int)(_track_ID[itrk])) {
              double thisClusterE = (_clusters_calo[fClusterAlgo][icalo].at(iclus)).cluster_E;
              E_clusters_match[icalo]+=thisClusterE;
              if(thisClusterE > maxClusterE) {
                maxClusterE = thisClusterE;
                Emax_cluster_M02 = (_clusters_calo[fClusterAlgo][icalo].at(iclus)).cluster_M02;
                Emax_cluster_M20 = (_clusters_calo[fClusterAlgo][icalo].at(iclus)).cluster_M20;
              }
            }
          } // track match loop
        } // is cluster matched
      } // cluster loop
    } // calo loop

    for(int imcp = 0; imcp < _nMCPart; imcp++) {		
      if(_mcpart_ID[imcp] == 1 && _mcpart_ID_parent[imcp] == 0 && _mcpart_PDG[imcp] == 11) {
        E_MC = _mcpart_E[imcp];
      }
    } // MC particle loop
    epidTree->Fill();
  } // track loop
}




void write_epid_tree() {

  epidFile->cd();
  epidTree->Write();
  epidFile->Close();	

}

void reset_variables() {

  p_track = -99.;
  eta_track = -99.;
  pT_track = -99.;
  theta_track = -99.;
  phi_track = -99.;
  for(int i = 0; i < maxcalo; i++) {
    E_clusters_all[i] = 0;
    E_clusters_match[i] = 0;
  }
  Emax_cluster_M02 = -99.;
  Emax_cluster_M20 = -99.;
  E_MC = -99.;
}
