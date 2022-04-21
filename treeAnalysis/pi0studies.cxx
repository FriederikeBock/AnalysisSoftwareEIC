#include <TLorentzVector.h>
// ANCHOR debug output verbosity
//Int_t verbosityERH = 0;

// ANCHOR create histograms globally
TH2F*  h_Etrue_Minv[_active_calo][_active_algo];
TH2F*  h_Erec_Minv[_active_calo][_active_algo];
TH2F*  h_Etrue_Minv_allcls[_active_calo][_active_algo];
TH2F*  h_Erec_Minv_allcls[_active_calo][_active_algo];
TH2F*  h_Etrue_Minv_allcalo;
TH2F*  h_Erec_Minv_allcalo;
TH2F*  h_Etrue_Minv_allcalo_unmatched;
TH2F*  h_Erec_Minv_allcalo_unmatched;
TH1F*  h_NClus[_active_calo][_active_algo];
TH1F*  h_Pi0Spec_total[_active_calo][_active_algo];
TH1F*  h_Pi0Spec_separatecls[_active_calo][_active_algo];
TH1F*  h_Pi0Spec_mergedclus[_active_calo][_active_algo];
TH2F*  h_Pi0Spec_total_Eta[_active_calo][_active_algo];
TH2F*  h_Pi0Spec_separatecls_Eta[_active_calo][_active_algo];
TH2F*  h_Pi0Spec_mergedclus_Eta[_active_calo][_active_algo];


void pi0studies(){
  std::vector<clustersStrct> _clusters_allcalo;

  CheckBarrelCaloVersion();

  for(int icalo=0;icalo<_active_calo;icalo++){

    //=====Loop only over enabled calorimeters ===========//
    if (!caloEnabled[icalo]) continue;

    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(!loadClusterizerInput(ialgo,icalo )) continue;

      int nbinsx = 500;
      // create output histograms
      if(!h_Etrue_Minv_allcalo) h_Etrue_Minv_allcalo  = new TH2F(Form("h_Etrue_Minv_allcalo"), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Erec_Minv_allcalo) h_Erec_Minv_allcalo  = new TH2F(Form("h_Erec_Minv_allcalo"), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Etrue_Minv_allcalo_unmatched) h_Etrue_Minv_allcalo_unmatched  = new TH2F(Form("h_Etrue_Minv_allcalo_unmatched"), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Erec_Minv_allcalo_unmatched) h_Erec_Minv_allcalo_unmatched  = new TH2F(Form("h_Erec_Minv_allcalo_unmatched"), "", nbinsx, 0, 50, nbinsx, 0, 0.5);

      if(!h_Etrue_Minv_allcls[icalo][ialgo]) h_Etrue_Minv_allcls[icalo][ialgo]  = new TH2F(Form("h_Etrue_Minv_allcls_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Erec_Minv_allcls[icalo][ialgo]) h_Erec_Minv_allcls[icalo][ialgo]  = new TH2F(Form("h_Erec_Minv_allcls_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Etrue_Minv[icalo][ialgo]) h_Etrue_Minv[icalo][ialgo]  = new TH2F(Form("h_Etrue_Minv_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_Erec_Minv[icalo][ialgo]) h_Erec_Minv[icalo][ialgo]  = new TH2F(Form("h_Erec_Minv_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", nbinsx, 0, 50, nbinsx, 0, 0.5);
      if(!h_NClus[icalo][ialgo]) h_NClus[icalo][ialgo]  = new TH1F(Form("h_NClus_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 50, 0, 50);
      if(!h_Pi0Spec_total[icalo][ialgo]) h_Pi0Spec_total[icalo][ialgo]  = new TH1F(Form("h_Pi0Spec_total_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 200, 0, 50);
      if(!h_Pi0Spec_separatecls[icalo][ialgo]) h_Pi0Spec_separatecls[icalo][ialgo]  = new TH1F(Form("h_Pi0Spec_separatecls_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 200, 0, 50);
      if(!h_Pi0Spec_mergedclus[icalo][ialgo]) h_Pi0Spec_mergedclus[icalo][ialgo]  = new TH1F(Form("h_Pi0Spec_mergedclus_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 200, 0, 50);
      if(!h_Pi0Spec_total_Eta[icalo][ialgo]) h_Pi0Spec_total_Eta[icalo][ialgo]  = new TH2F(Form("h_Pi0Spec_total_Eta_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 200, 0, 50, 160, -4, 4);
      if(!h_Pi0Spec_separatecls_Eta[icalo][ialgo]) h_Pi0Spec_separatecls_Eta[icalo][ialgo]  = new TH2F(Form("h_Pi0Spec_separatecls_Eta_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 200, 0, 50, 160, -4, 4);
      if(!h_Pi0Spec_mergedclus_Eta[icalo][ialgo]) h_Pi0Spec_mergedclus_Eta[icalo][ialgo]  = new TH2F(Form("h_Pi0Spec_mergedclus_Eta_%s_%s", str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 200, 0, 50, 160, -4, 4);

      // calculations only necessary if at least two clusters are present
      h_NClus[icalo][ialgo]->Fill(_clusters_calo[ialgo][icalo].size());
      // if(_clusters_calo[ialgo][icalo].size() > 1){
      //   bool processclusters = true;
      //   // if(icalo==kFEMC){
      //   //   if((_clusters_calo[ialgo][icalo].at(0)).cluster_Eta < 2.0) processclusters = false;
      //   // } else if(icalo==kBECAL){
      //   //   if(abs((_clusters_calo[ialgo][icalo].at(0)).cluster_Eta) > 0.5) processclusters = false;
      //   // } else if(icalo==kEEMCG){
      //   //   if(abs((_clusters_calo[ialgo][icalo].at(0)).cluster_Eta) > -2.0) processclusters = false;
      //   // }
      //   if(processclusters){

      //     // try with x,y,z of cluster relative to smeared vertex
      //     // float _vertex_x;
      //     // float _vertex_y;
      //     // float _vertex_z;
      //     TVector3 vecClus1((_clusters_calo[ialgo][icalo].at(0)).cluster_X - _vertex_x,(_clusters_calo[ialgo][icalo].at(0)).cluster_Y - _vertex_y,(_clusters_calo[ialgo][icalo].at(0)).cluster_Z - _vertex_z);
      //     // cout << vecClus1.Eta() << "\t" << (_clusters_calo[ialgo][icalo].at(0)).cluster_Eta << "\t" << _vertex_z << endl;

      //     // try hybrid calo methods

      //     //pi0 reco efficiency

      //     double Erec = (_clusters_calo[ialgo][icalo].at(0)).cluster_E + (_clusters_calo[ialgo][icalo].at(1)).cluster_E ;
      //     double E1 = (_clusters_calo[ialgo][icalo].at(0)).cluster_E;
      //     double E2 = (_clusters_calo[ialgo][icalo].at(1)).cluster_E;
      //     double theta1 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(0)).cluster_Eta));
      //     double theta2 = 2*atan(exp(-(_clusters_calo[ialgo][icalo].at(1)).cluster_Eta));
      //     double phi1   = (_clusters_calo[ialgo][icalo].at(0)).cluster_Phi;
      //     double phi2   = (_clusters_calo[ialgo][icalo].at(1)).cluster_Phi;

      //     double pX1    = E1*sin(theta1)*cos(phi1);
      //     double pY1    = E1*sin(theta1)*sin(phi1);
      //     double pZ1    = E1*cos(theta1);
      //     TLorentzVector photon1(pX1, pY1, pZ1, E1);

      //     double pX2    = E2*sin(theta2)*cos(phi2);
      //     double pY2    = E2*sin(theta2)*sin(phi2);
      //     double pZ2    = E2*cos(theta2);
      //     TLorentzVector photon2(pX2, pY2, pZ2, E2);

      //     TLorentzVector pi0 = photon1 + photon2;
      //     double Minv = pi0.M();
      //     //cout << "====== M ========" << Minv << " " << photon1.M() << " " << photon1.Eta() << " " <<  (_clusters_calo[ialgo][icalo].at(0)).cluster_Eta  << " " << photon1.Phi() << " " << phi1 << endl;
      //     //cout << "====== M ========" << Minv << " " << photon2.M() << " " << photon2.Eta() << " " <<  (_clusters_calo[ialgo][icalo].at(1)).cluster_Eta  << " " << photon2.Phi() << " " << phi2 << endl;

      //     h_Erec_Minv[icalo][ialgo]->Fill(Erec,Minv);
      //     // if(_mcpart_PDG[0]!=111) cout << "not a pion!" << endl;
      //     // cout << _mcpart_E[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID] << endl;
      //     if(hepMCavailable){
      //       int bcid_clus1 = _mcpart_BCID[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID]-1; // -1 due to array filling
      //       int bcid_clus2 = _mcpart_BCID[(_clusters_calo[ialgo][icalo].at(1)).cluster_trueID]-1; // -1 due to array filling
      //       if(_hepmcp_m1[bcid_clus1]!=0){
      //           // cout << "cls1\thepPDG \tMCPDG \tMotherPDG " << endl;
      //           // cout << "\t" << _hepmcp_PDG[bcid_clus1] << "\t " << _mcpart_PDG[(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_trueID]<< "\t " << _hepmcp_PDG[_hepmcp_m1[bcid_clus1]-1] << endl;
      //         if(_hepmcp_PDG[_hepmcp_m1[bcid_clus1]-1]==111){
      //           h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_hepmcp_E[_hepmcp_m1[bcid_clus1]-1],Minv);
      //         }
      //       }
      //       if(_hepmcp_m1[bcid_clus2]!=0){
      //           // cout << "cls2\thepPDG \tMCPDG \tMotherPDG " << endl;
      //           // cout << "\t" << _hepmcp_PDG[bcid_clus2] << "\t " << _mcpart_PDG[(_clusters_calo[ialgo][icalo].at(iclus2)).cluster_trueID]<< "\t " << _hepmcp_PDG[_hepmcp_m1[bcid_clus2]-1] << endl;
      //         if(_hepmcp_PDG[_hepmcp_m1[bcid_clus2]-1]==111){
      //           h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_hepmcp_E[_hepmcp_m1[bcid_clus2]-1],Minv);
      //         }
      //       }
      //     } else {
      //       if(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID]==111){
      //         h_Etrue_Minv[icalo][ialgo]->Fill(_mcpart_E[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID],Minv);
      //       } else if(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(1)).cluster_trueID]==111){
      //         h_Etrue_Minv[icalo][ialgo]->Fill(_mcpart_E[(_clusters_calo[ialgo][icalo].at(1)).cluster_trueID],Minv);
      //       }
      //     }
      //     // get true pi0 candidates, as well as candidates with only one true contribution
      //     bool passedMCcheck = false;
      //     if(passedMCcheck){

      //     }
      //   }
      // }

      for(Int_t iclus1=0; iclus1<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus1++){
        if(ialgo==kMA)_clusters_allcalo.push_back(_clusters_calo[ialgo][icalo].at(iclus1));
        // if((_clusters_calo[ialgo][icalo].at(iclus1)).cluster_NTowers<2) continue;
        // if(icalo==kFEMC){
        //   if((_clusters_calo[ialgo][icalo].at(iclus1)).cluster_Eta < 2.5) continue;
        // } else if(icalo==kBECAL){
        //   if(abs((_clusters_calo[ialgo][icalo].at(iclus1)).cluster_Eta) > 0.8) continue;
        // } else if(icalo==kEEMCG){
        //   if((_clusters_calo[ialgo][icalo].at(iclus1)).cluster_Eta > -2.5) continue;
        // }
        for(Int_t iclus2=iclus1+1; iclus2<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus2++){
          // if(icalo==kFEMC){
          //   if((_clusters_calo[ialgo][icalo].at(iclus2)).cluster_Eta < 2.5) continue;
          // } else if(icalo==kBECAL){
          //   if(abs((_clusters_calo[ialgo][icalo].at(iclus2)).cluster_Eta) > 0.8) continue;
          // } else if(icalo==kEEMCG){
          //   if(abs((_clusters_calo[ialgo][icalo].at(iclus2)).cluster_Eta) > -2.5) continue;
          // }
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
          // h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_mcpart_E[(_clusters_calo[ialgo][icalo].at(0)).cluster_trueID],Minv);
          if(hepMCavailable){
            int bcid_clus1 = _mcpart_BCID[(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_trueID]-1; // -1 due to array filling
            int bcid_clus2 = _mcpart_BCID[(_clusters_calo[ialgo][icalo].at(iclus2)).cluster_trueID]-1; // -1 due to array filling
            if(_hepmcp_m1[bcid_clus1]!=0){
                // cout << "cls1\thepPDG \tMCPDG \tMotherPDG " << endl;
                // cout << "\t" << _hepmcp_PDG[bcid_clus1] << "\t " << _mcpart_PDG[(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_trueID]<< "\t " << _hepmcp_PDG[_hepmcp_m1[bcid_clus1]-1] << endl;
              if(_hepmcp_PDG[_hepmcp_m1[bcid_clus1]-1]==111){
                h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_hepmcp_E[_hepmcp_m1[bcid_clus1]-1],Minv);
              }
            }
            if(_hepmcp_m1[bcid_clus2]!=0){
                // cout << "cls2\thepPDG \tMCPDG \tMotherPDG " << endl;
                // cout << "\t" << _hepmcp_PDG[bcid_clus2] << "\t " << _mcpart_PDG[(_clusters_calo[ialgo][icalo].at(iclus2)).cluster_trueID]<< "\t " << _hepmcp_PDG[_hepmcp_m1[bcid_clus2]-1] << endl;
              if(_hepmcp_PDG[_hepmcp_m1[bcid_clus2]-1]==111){
                h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_hepmcp_E[_hepmcp_m1[bcid_clus2]-1],Minv);
              }
            }
          } else {
            // cout << _mcpart_PDG[(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_trueID] << "\tmother: " << _mcpart_PDG[_mcpart_ID_parent[(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_trueID]] << endl;
            if(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_trueID]==111){
              h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_mcpart_E[(_clusters_calo[ialgo][icalo].at(iclus1)).cluster_trueID],Minv);
            } else if(_mcpart_PDG[(_clusters_calo[ialgo][icalo].at(iclus2)).cluster_trueID]==111){
              h_Etrue_Minv_allcls[icalo][ialgo]->Fill(_mcpart_E[(_clusters_calo[ialgo][icalo].at(iclus2)).cluster_trueID],Minv);
            }
          }
        } // iclus2 loop
      } // iclus1 loop


      if(!hepMCavailable){ // this means it is a single particle simulation
      // cout << "new evt" << endl;
        for(int imc=0; imc<_nMCPart; imc++){
          // approach is to find the cluster(s) for each generated pi0
          // cout << _mcpart_PDG[imc] << "motherID: " << _mcpart_ID_parent[imc]<< endl;
          if(_mcpart_PDG[imc]==111){
            // TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
            float trueetaPi0 = _mcpart_Eta[imc];
            if((trueetaPi0<-2.5 || trueetaPi0>-1.5) && icalo==kEEMC) continue;
            if((trueetaPi0<-1.2 || trueetaPi0>1.2) && icalo==kBECAL) continue;
            if((trueetaPi0<1.6 || trueetaPi0>2.5) && icalo==kFEMC) continue;
            if((trueetaPi0<2.0 || trueetaPi0>3.0) && icalo==kDRCALO) continue;

            int gammasInAcc = 0;
            TVector3 g1,g2;
            for(int imc2=0; imc2<_nMCPart; imc2++){
              if(_mcpart_ID_parent[imc2]==1 && _mcpart_PDG[imc2]==22){
                // TVector3 mcpartvec2(_mcpart_px[imc2],_mcpart_py[imc2],_mcpart_pz[imc2]);
                float trueetaGamma = _mcpart_Eta[imc2];
                if((trueetaGamma<-2.5 || trueetaGamma>-1.5) && icalo==kEEMC) continue;
                if((trueetaGamma<-1.2 || trueetaGamma>1.2) && icalo==kBECAL) continue;
                if((trueetaGamma<1.6 || trueetaGamma>2.5) && icalo==kFEMC) continue;
                if((trueetaGamma<2.0 || trueetaGamma>3.0) && icalo==kDRCALO) continue;
                gammasInAcc++;
                if(gammasInAcc==1) g1.SetXYZ(_mcpart_px[imc2],_mcpart_py[imc2],_mcpart_pz[imc2]);
                else if(gammasInAcc==2) g2.SetXYZ(_mcpart_px[imc2],_mcpart_py[imc2],_mcpart_pz[imc2]);
              }
            }
            // both decay photons must be in calorimeter acceptance
            if(gammasInAcc!=2){
              // cout << "only " << gammasInAcc << " found" << endl;
              continue;
            }
            bool ismerged = false;
            float caloPos = 0;
            switch (icalo){
              case kFEMC:
                caloPos = 315;
                g1 *= caloPos/g1.Z();
                g2 *= caloPos/g2.Z();
                if((g1-g2).Mag()<sqrt(2*pow(2.5,2))) ismerged = true;
                break;
              case kEEMC:
                caloPos = -240;
                g1 *= caloPos/g1.Z();
                g2 *= caloPos/g2.Z();
                if((g1-g2).Mag()<sqrt(2*pow(3.2,2))) ismerged = true;
                break;
              case kBECAL:
                caloPos = 64;
                g1 *= caloPos/sqrt(pow(g1.X(),2)+pow(g1.Y(),2));
                g2 *= caloPos/sqrt(pow(g2.X(),2)+pow(g2.Y(),2));
                if((g1-g2).Mag()<sqrt(2*pow(4.0,2))) ismerged = true;
                break;
              case kDRCALO:
                caloPos = 270;
                g1 *= caloPos/g1.Z();
                g2 *= caloPos/g2.Z();
                if((g1-g2).Mag()<1.0) ismerged = true;
                break;
              default:
                ismerged = false;
                break;
            }

            if(ismerged){
              h_Pi0Spec_total[ialgo][icalo]->Fill(_mcpart_E[imc]);
              h_Pi0Spec_mergedclus[ialgo][icalo]->Fill(_mcpart_E[imc]);
              h_Pi0Spec_total_Eta[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
              h_Pi0Spec_mergedclus_Eta[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
            } else {
              h_Pi0Spec_total[ialgo][icalo]->Fill(_mcpart_E[imc]);
              h_Pi0Spec_separatecls[ialgo][icalo]->Fill(_mcpart_E[imc]);
              h_Pi0Spec_total_Eta[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
              h_Pi0Spec_separatecls[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
            }
            // int iclusfound = 0;
            // for(Int_t iclus1=0; iclus1<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus1++){
            //   if((_clusters_calo[ialgo][icalo].at(iclus1)).cluster_E > (0.05*_mcpart_E[imc])){
            //     iclusfound++;
            //   }
            // }
            // if(iclusfound==1){
            //   h_Pi0Spec_total[ialgo][icalo]->Fill(_mcpart_E[imc]);
            //   h_Pi0Spec_mergedclus[ialgo][icalo]->Fill(_mcpart_E[imc]);
            //   h_Pi0Spec_total_Eta[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
            //   h_Pi0Spec_mergedclus_Eta[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
            // } else if(iclusfound>1){
            //   h_Pi0Spec_total[ialgo][icalo]->Fill(_mcpart_E[imc]);
            //   h_Pi0Spec_separatecls[ialgo][icalo]->Fill(_mcpart_E[imc]);
            //   h_Pi0Spec_total_Eta[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
            //   h_Pi0Spec_separatecls[ialgo][icalo]->Fill(_mcpart_E[imc],trueetaPi0);
            // }
          }
        }
      }

    } // ialgo loop
  } // icalo loop

  // hybrid method with all calos
  for(Int_t iclus1=0; iclus1<(Int_t)_clusters_allcalo.size(); iclus1++){
    // if((_clusters_allcalo.at(iclus1)).cluster_NTowers<2) continue;
    for(Int_t iclus2=iclus1+1; iclus2<(Int_t)_clusters_allcalo.size(); iclus2++){
      // if((_clusters_allcalo.at(iclus2)).cluster_NTowers<2) continue;
      double Erec = (_clusters_allcalo.at(iclus1)).cluster_E + (_clusters_allcalo.at(iclus2)).cluster_E ;
      double E1 = (_clusters_allcalo.at(iclus1)).cluster_E;
      double E2 = (_clusters_allcalo.at(iclus2)).cluster_E;
      double theta1 = 2*atan(exp(-(_clusters_allcalo.at(iclus1)).cluster_Eta));
      double theta2 = 2*atan(exp(-(_clusters_allcalo.at(iclus2)).cluster_Eta));
      double phi1   = (_clusters_allcalo.at(iclus1)).cluster_Phi;
      double phi2   = (_clusters_allcalo.at(iclus2)).cluster_Phi;

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

      h_Erec_Minv_allcalo->Fill(Erec,Minv);

      if(hepMCavailable){
        int bcid_clus1 = _mcpart_BCID[(_clusters_allcalo.at(iclus1)).cluster_trueID]-1; // -1 due to array filling
        int bcid_clus2 = _mcpart_BCID[(_clusters_allcalo.at(iclus2)).cluster_trueID]-1; // -1 due to array filling
        if(_hepmcp_m1[bcid_clus1]!=0){
            // cout << "cls1\thepPDG \tMCPDG \tMotherPDG " << endl;
            // cout << "\t" << _hepmcp_PDG[bcid_clus1] << "\t " << _mcpart_PDG[(_clusters_allcalo.at(iclus1)).cluster_trueID]<< "\t " << _hepmcp_PDG[_hepmcp_m1[bcid_clus1]-1] << endl;
          if(_hepmcp_PDG[_hepmcp_m1[bcid_clus1]-1]==111){
            h_Etrue_Minv_allcalo->Fill(_hepmcp_E[_hepmcp_m1[bcid_clus1]-1],Minv);
          }
        }
        if(_hepmcp_m1[bcid_clus2]!=0){
            // cout << "cls2\thepPDG \tMCPDG \tMotherPDG " << endl;
            // cout << "\t" << _hepmcp_PDG[bcid_clus2] << "\t " << _mcpart_PDG[(_clusters_allcalo.at(iclus2)).cluster_trueID]<< "\t " << _hepmcp_PDG[_hepmcp_m1[bcid_clus2]-1] << endl;
          if(_hepmcp_PDG[_hepmcp_m1[bcid_clus2]-1]==111){
            h_Etrue_Minv_allcalo->Fill(_hepmcp_E[_hepmcp_m1[bcid_clus2]-1],Minv);
          }
        }
      } else {
        if(_mcpart_PDG[(_clusters_allcalo.at(iclus1)).cluster_trueID]==111){
          h_Etrue_Minv_allcalo->Fill(_mcpart_E[(_clusters_allcalo.at(iclus1)).cluster_trueID],Minv);
        } else if(_mcpart_PDG[(_clusters_allcalo.at(iclus2)).cluster_trueID]==111){
          h_Etrue_Minv_allcalo->Fill(_mcpart_E[(_clusters_allcalo.at(iclus2)).cluster_trueID],Minv);
        }
      }
    }
  }

  // hybrid method with all calos only neutral clusters!
  for(Int_t iclus1=0; iclus1<(Int_t)_clusters_allcalo.size(); iclus1++){
    // if((_clusters_allcalo.at(iclus1)).cluster_NTowers<2) continue;
    if((_clusters_allcalo.at(iclus1)).cluster_isMatched) continue;
    for(Int_t iclus2=iclus1+1; iclus2<(Int_t)_clusters_allcalo.size(); iclus2++){
    if((_clusters_allcalo.at(iclus2)).cluster_isMatched) continue;
      // if((_clusters_allcalo.at(iclus2)).cluster_NTowers<2) continue;
      double Erec = (_clusters_allcalo.at(iclus1)).cluster_E + (_clusters_allcalo.at(iclus2)).cluster_E ;
      double E1 = (_clusters_allcalo.at(iclus1)).cluster_E;
      double E2 = (_clusters_allcalo.at(iclus2)).cluster_E;
      double theta1 = 2*atan(exp(-(_clusters_allcalo.at(iclus1)).cluster_Eta));
      double theta2 = 2*atan(exp(-(_clusters_allcalo.at(iclus2)).cluster_Eta));
      double phi1   = (_clusters_allcalo.at(iclus1)).cluster_Phi;
      double phi2   = (_clusters_allcalo.at(iclus2)).cluster_Phi;

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

      h_Erec_Minv_allcalo_unmatched->Fill(Erec,Minv);
      if(_mcpart_PDG[(_clusters_allcalo.at(iclus1)).cluster_trueID]==111){
        h_Etrue_Minv_allcalo_unmatched->Fill(_mcpart_E[(_clusters_allcalo.at(iclus1)).cluster_trueID],Minv);
      } else if(_mcpart_PDG[(_clusters_allcalo.at(iclus2)).cluster_trueID]==111){
        h_Etrue_Minv_allcalo_unmatched->Fill(_mcpart_E[(_clusters_allcalo.at(iclus2)).cluster_trueID],Minv);
      }
    }
  }
}

// ANCHOR save function after event loop

void pi0studiesSave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_pi0.root",outputDir.Data()),"RECREATE");

  // write histograms
  if(h_Etrue_Minv_allcalo)  h_Etrue_Minv_allcalo->Write();
  if(h_Erec_Minv_allcalo) h_Erec_Minv_allcalo->Write();
  if(h_Etrue_Minv_allcalo_unmatched)  h_Etrue_Minv_allcalo_unmatched->Write();
  if(h_Erec_Minv_allcalo_unmatched) h_Erec_Minv_allcalo_unmatched->Write();
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
      if(h_Pi0Spec_total[icalo][ialgo]) h_Pi0Spec_total[icalo][ialgo]->Write();
      if(h_Pi0Spec_separatecls[icalo][ialgo]) h_Pi0Spec_separatecls[icalo][ialgo]->Write();
      if(h_Pi0Spec_mergedclus[icalo][ialgo]) h_Pi0Spec_mergedclus[icalo][ialgo]->Write();
      if(h_Pi0Spec_total_Eta[icalo][ialgo]) h_Pi0Spec_total_Eta[icalo][ialgo]->Write();
      if(h_Pi0Spec_separatecls_Eta[icalo][ialgo]) h_Pi0Spec_separatecls_Eta[icalo][ialgo]->Write();
      if(h_Pi0Spec_mergedclus_Eta[icalo][ialgo]) h_Pi0Spec_mergedclus_Eta[icalo][ialgo]->Write();
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
