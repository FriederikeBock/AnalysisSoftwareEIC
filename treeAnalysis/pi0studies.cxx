
// ANCHOR debug output verbosity
//Int_t verbosityERH = 0;

// ANCHOR create histograms globally
TH2F*  h_Etrue_Minv[_active_calo][_active_algo];
TH2F*  h_Erec_Minv[_active_calo][_active_algo];
TH2F*  h_Etrue_Minv_allcls[_active_calo][_active_algo];
TH2F*  h_Erec_Minv_allcls[_active_calo][_active_algo];
TH1F*  h_NClus[_active_calo][_active_algo];

void pi0studies(){

  CheckBarrelCaloVersion();

  for(int icalo=0;icalo<_active_calo;icalo++){

    //=====Loop only over enabled calorimeters ===========//
    if (!caloEnabled[icalo]) continue;

    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!loadClusterizerInput(ialgo,icalo )) continue;

      std::sort(_clusters_calo[ialgo][icalo].begin(), _clusters_calo[ialgo][icalo].end(), &acompareCl);

      int nbinsx = 500;
      // create output histograms
      if(!h_Etrue_Minv_allcls[icalo][ialgo]) h_Etrue_Minv_allcls[icalo][ialgo]  = new TH2F(Form("h_Etrue_Minv_allcls_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Erec_Minv_allcls[icalo][ialgo]) h_Erec_Minv_allcls[icalo][ialgo]  = new TH2F(Form("h_Erec_Minv_allcls_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Etrue_Minv[icalo][ialgo]) h_Etrue_Minv[icalo][ialgo]  = new TH2F(Form("h_Etrue_Minv_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Erec_Minv[icalo][ialgo]) h_Erec_Minv[icalo][ialgo]  = new TH2F(Form("h_Erec_Minv_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_NClus[icalo][ialgo]) h_NClus[icalo][ialgo]  = new TH1F(Form("h_NClus_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 50, 0, 50);

      // calculations only necessary if at least two clusters are present
      h_NClus[icalo][ialgo]->Fill(_clusters_calo[ialgo][icalo].size());
      if(_clusters_calo[ialgo][icalo].size() > 1){
        bool processclusters = true;
        // if(icalo==kFEMC){
        //   if((_clusters_calo[ialgo][icalo].at(0)).cluster_Eta < 2.0) processclusters = false;
        // } else if(icalo==kBECAL){
        //   if(abs((_clusters_calo[ialgo][icalo].at(0)).cluster_Eta) > 0.5) processclusters = false;
        // } else if(icalo==kEEMCG){
        //   if(abs((_clusters_calo[ialgo][icalo].at(0)).cluster_Eta) > -2.0) processclusters = false;
        // }
        if(processclusters){

          // try with x,y,z of cluster relative to smeared vertex

          // try hybrid calo methods

          //pi0 reco efficiency

          double Erec = (_clusters_calo[ialgo][icalo].at(0)).cluster_E + (_clusters_calo[ialgo][icalo].at(1)).cluster_E ;
          double E1 = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
          double E2 = (_clusters_calo[ialgo][icalo].at(1)).cluster_E;
          double theta1 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(0)).cluster_Eta));
          double theta2 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(1)).cluster_Eta));
          double phi1   = (_clusters_calo[ialgo][icalo].at(0)).cluster_Phi;
          double phi2   = (_clusters_calo[ialgo][icalo].at(1)).cluster_Phi;

          double pX1    = E1*sin(theta1)*cos(phi1);
          double pY1    = E1*sin(theta1)*sin(phi1);
          double pZ1    = E1*cos(theta1);
          TLorentzVector photon1(pX1, pY1, pZ1, E1);

          double pX2    = E2*sin(theta2)*cos(phi2);
          double pY2    = E2*sin(theta2)*sin(phi2);
          double pZ2    = E2*cos(theta2);
          TLorentzVector photon2(pX2, pY2, pZ2, E2);

          TLorentzVector pi0 = photon1 + photon2;
          double Minv = pi0.M();
          //cout << "====== M ========" << Minv << " " << photon1.M() << " " << photon1.Eta() << " " <<  (_clusters_calo[ialgo][icalo].at(0)).cluster_Eta  << " " << photon1.Phi() << " " << phi1 << endl;
          //cout << "====== M ========" << Minv << " " << photon2.M() << " " << photon2.Eta() << " " <<  (_clusters_calo[ialgo][icalo].at(1)).cluster_Eta  << " " << photon2.Phi() << " " << phi2 << endl;

          h_Erec_Minv[icalo][ialgo]->Fill(Erec,Minv);
          // if(_mcpart_PDG[0]!=111) cout << "not a pion!" << endl;
          // cout << _mcpart_E[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID] << endl;
          h_Etrue_Minv[icalo][ialgo]->Fill(_mcpart_E[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID],Minv);

          // get true pi0 candidates, as well as candidates with only one true contribution
          bool passedMCcheck = false;
          if(passedMCcheck){

          }
        }
      }

      for(Int_t iclus1=0; iclus1<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus1++){
        // if((_clusters_calo[ialgo][icalo].at(iclus1)).cluster_NTowers<2) continue;
        for(Int_t iclus2=iclus1+1; iclus2<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus2++){
          // if((_clusters_calo[ialgo][icalo].at(iclus2)).cluster_NTowers<2) continue;
          double Erec = (_clusters_calo[ialgo][icalo].at(iclus1)).cluster_E + (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_E ;
          double E1 = (_clusters_calo[ialgo][icalo].at(iclus1)).cluster_E;
          double E2 = (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_E;
          double theta1 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_Eta));
          double theta2 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(iclus2)).cluster_Eta));
          double phi1   = (_clusters_calo[ialgo][icalo].at(iclus1)).cluster_Phi;
          double phi2   = (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_Phi;

          double pX1    = E1*sin(theta1)*cos(phi1);
          double pY1    = E1*sin(theta1)*sin(phi1);
          double pZ1    = E1*cos(theta1);
          TLorentzVector photon1(pX1, pY1, pZ1, E1);

          double pX2    = E2*sin(theta2)*cos(phi2);
          double pY2    = E2*sin(theta2)*sin(phi2);
          double pZ2    = E2*cos(theta2);
          TLorentzVector photon2(pX2, pY2, pZ2, E2);

          TLorentzVector pi0 = photon1 + photon2;
          double Minv = pi0.M();

          h_Erec_Minv_allcls[icalo][ialgo]->Fill(Erec,Minv);
          h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_mcpart_E[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID],Minv);
        }
      }
    }
  }
}

// ANCHOR save function after event loop

void pi0studiesSave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_pi0.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      //====== Eta =====================================================/
      if(h_NClus[icalo][ialgo])  h_NClus[icalo][ialgo]->Write();
      if(h_Etrue_Minv[icalo][ialgo])  h_Etrue_Minv[icalo][ialgo]->Write();
      if(h_Erec_Minv[icalo][ialgo]) h_Erec_Minv[icalo][ialgo]->Write();
      if(h_Etrue_Minv_allcls[icalo][ialgo])  h_Etrue_Minv_allcls[icalo][ialgo]->Write();
      if(h_Erec_Minv_allcls[icalo][ialgo]) h_Erec_Minv_allcls[icalo][ialgo]->Write();
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
