
// ANCHOR debug output verbosity
Int_t verbosityRH = 0;

const int particleSpRes = 7;
TString addNameRes[particleSpRes]   = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_" ,""};
int pdgCodeRes[particleSpRes]       = {22, 11, 211, 2212, 321, 2112, -1};

// ANCHOR create histograms globally
TH2F*  h_RH_Reso_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_Reso_smear_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_Reso_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoCalib_smeared_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_Reso_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];

TH2F*  h_RH_Reso_Ehighest[particleSpRes][_active_calo][_active_algo];


TH2F*  h_RH_ResoComb_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_smear_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoCombCalib_smeared_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_RH_ResoComb_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];

// ANCHOR main function to be called in event loop
void resolutionhistos(){

  // f_FEMC_calib->SetParameters(paramsFEMCcalib);

  // f_FHCAL_calib->SetParameters(paramsHCALcalib);

  CheckBarrelCaloVersion();
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(verbosityRH) cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << endl;
      if(!loadClusterizerInput( ialgo, icalo )) continue;
      
      std::sort(_clusters_calo[ialgo][icalo].begin(), _clusters_calo[ialgo][icalo].end(), &acompareCl);

      
      int nbinsx = 400;
      if(icalo==kFHCAL) nbinsx = 200;
      if(icalo==kEHCAL || icalo==kHCALOUT || icalo==kHCALOUT) nbinsx = 100;
      if(icalo==kLFHCAL) nbinsx = 200;
      if(icalo==kEEMC) nbinsx = 200;
      if(icalo==kBECAL) nbinsx = 200;
      if(icalo==kFEMC) nbinsx = 200;
      
      for(Int_t sp = 0; sp< particleSpRes; sp++){
        if(!h_RH_Reso_E[sp][icalo][ialgo]) h_RH_Reso_E[sp][icalo][ialgo]              = new TH2F(Form("h_RH_Reso_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        if(!h_RH_Reso_Ehighest[sp][icalo][ialgo]) h_RH_Reso_Ehighest[sp][icalo][ialgo] = new TH2F(Form("h_RH_Reso_%sEhighest_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        if(!h_RH_Reso_smear_E[sp][icalo][ialgo]) h_RH_Reso_smear_E[sp][icalo][ialgo]  = new TH2F(Form("h_RH_Reso_%ssmear_E_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        if(!h_RH_Reso_Erec[sp][icalo][ialgo])h_RH_Reso_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_RH_Reso_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        if(!h_RH_ResoCalib_E[sp][icalo][ialgo])h_RH_ResoCalib_E[sp][icalo][ialgo]      = new TH2F(Form("h_RH_ResoCalib_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
      
        if(!h_RH_ResoCalib_Erec[sp][icalo][ialgo])h_RH_ResoCalib_Erec[sp][icalo][ialgo] = new TH2F(Form("h_RH_ResoCalib_%sErec_%s_%s", addNameRes[sp].Data(),  str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        if(!h_RH_ResoCalib_smeared_Erec[sp][icalo][ialgo])h_RH_ResoCalib_smeared_Erec[sp][icalo][ialgo] = new TH2F(Form("h_RH_ResoCalib_smeared_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
      
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(!h_RH_Reso_E_Eta[sp][icalo][ialgo][eT])h_RH_Reso_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_RH_Reso_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        }        
      }
      
      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus++){
        // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
        // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

        // cluster should have at least 2 towers
        if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers<2) continue;

        float clusterE_calo = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E;
        clusterE_calo/=getCalibrationValue((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, icalo, ialgo);
        
        int mcID = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;
        // loop over MC particles

        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[mcID])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
        
        // find true MC particle for given cluster
        // if(verbosityRH>1) cout << "\tfound MC:  particle " << mcID << "\twith E = " << _mcpart_E[mcID] << " GeV" << endl;
        h_RH_Reso_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);
        if (pClass != -1) h_RH_Reso_E[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);

        if (iclus == 0 && _nMCPart==1 ){
          h_RH_Reso_Ehighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);
          if (pClass != -1) h_RH_Reso_Ehighest[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);
        
        }
        float smearedClusE = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E*=getEnergySmearing(icalo,ialgo);
        h_RH_Reso_smear_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],smearedClusE/_mcpart_E[mcID]);
        if (pClass != -1) h_RH_Reso_smear_E[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],smearedClusE/_mcpart_E[mcID]);

        Int_t et = 0;
        while ( ( (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta > partEta[et+1] ) && ( et < nEta )) et++;
        h_RH_Reso_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);
        if (pClass != -1) h_RH_Reso_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);
        
        h_RH_Reso_Erec[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);
        if (pClass != -1) h_RH_Reso_Erec[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID]);
        
        h_RH_ResoCalib_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],clusterE_calo/_mcpart_E[mcID]);
        if (pClass != -1) h_RH_ResoCalib_E[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],clusterE_calo/_mcpart_E[mcID]);

        h_RH_ResoCalib_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[mcID]);
        if (pClass != -1) h_RH_ResoCalib_Erec[pClass][icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[mcID]);

        h_RH_ResoCalib_smeared_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterE_calo,(clusterE_calo*getEnergySmearing(icalo,ialgo))/_mcpart_E[mcID]);
        if (pClass != -1) h_RH_ResoCalib_smeared_Erec[pClass][icalo][ialgo]->Fill(clusterE_calo,(clusterE_calo*getEnergySmearing(icalo,ialgo))/_mcpart_E[mcID]);
      }
    }
  }
  
  if (_nMCPart==1){
    for(int icalo=0;icalo<_active_calo;icalo++){
      if (!caloEnabled[icalo]) continue;          // calorimeter is enabled
      if (_combCalo[icalo] == -1 ) continue;      // inner Hcal or Ecal defined
      if (!caloEnabled[_combCalo[icalo]]) continue; // inner HCal or Ecal ran
      for(int ialgo=0;ialgo<_active_algo;ialgo++){
        if (_combCalo[icalo] == -1 ) continue;
        bool bRefCalo = loadClusterizerInput( ialgo, icalo);
        bool b1stCalo = loadClusterizerInput( ialgo, _combCalo[icalo]);
        
        std::sort(_clusters_calo[ialgo][_combCalo[icalo]].begin(), _clusters_calo[ialgo][_combCalo[icalo]].end(), &acompareCl);    
          
        bool b2ndCalo = false;
        if (_combCalo2[icalo] != -1 ){
          if (caloEnabled[_combCalo2[icalo]]){ 
            b2ndCalo = loadClusterizerInput( ialgo, _combCalo2[icalo]);
          }
          std::sort(_clusters_calo[ialgo][_combCalo2[icalo]].begin(), _clusters_calo[ialgo][_combCalo2[icalo]].end(), &acompareCl);  
        }
        
        if (!(b1stCalo || b2ndCalo)) continue;
          
        if (verbosityRH){
          cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo ;
          if (b1stCalo)
            cout  << " with ECal: "<< _combCalo[icalo] ;
          if (b2ndCalo)
            cout << " and 2nd calo: " <<  _combCalo2[icalo] << endl;
          else 
            cout << endl;
        }
        int nbinsx = 400;
        if(icalo==kFHCAL) nbinsx = 200;
        if(icalo==kEHCAL || icalo==kHCALOUT || icalo==kHCALOUT) nbinsx = 100;
        if(icalo==kLFHCAL ) nbinsx = 200;
        if(icalo==kEEMC) nbinsx = 200;
        if(icalo==kBECAL) nbinsx = 200;
        if(icalo==kFEMC) nbinsx = 200;
        
        for(Int_t sp = 0; sp< particleSpRes; sp++){
          if(!h_RH_ResoComb_E[sp][icalo][ialgo]) h_RH_ResoComb_E[sp][icalo][ialgo]              = new TH2F(Form("h_RH_ResoComb_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
          if(!h_RH_ResoComb_smear_E[sp][icalo][ialgo]) h_RH_ResoComb_smear_E[sp][icalo][ialgo]  = new TH2F(Form("h_RH_ResoComb_%ssmear_E_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
          if(!h_RH_ResoComb_Erec[sp][icalo][ialgo])h_RH_ResoComb_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_RH_ResoComb_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
          if(!h_RH_ResoCombCalib_E[sp][icalo][ialgo])h_RH_ResoCombCalib_E[sp][icalo][ialgo]      = new TH2F(Form("h_RH_ResoCombCalib_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        
          if(!h_RH_ResoCombCalib_Erec[sp][icalo][ialgo])h_RH_ResoCombCalib_Erec[sp][icalo][ialgo] = new TH2F(Form("h_RH_ResoCombCalib_%sErec_%s_%s", addNameRes[sp].Data(),  str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
          if(!h_RH_ResoCombCalib_smeared_Erec[sp][icalo][ialgo])h_RH_ResoCombCalib_smeared_Erec[sp][icalo][ialgo] = new TH2F(Form("h_RH_ResoCombCalib_smeared_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsx, 0, 2);
        
          for (Int_t eT = 0; eT < nEta+1; eT++){
            if(!h_RH_ResoComb_E_Eta[sp][icalo][ialgo][eT])h_RH_ResoComb_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_RH_ResoComb_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsx, 0, 2);
          }        
        }

        Int_t highestECl = -1;
        if (bRefCalo){
          highestECl = 0;
        }
        Int_t highestECl2 = -1;
        if (b1stCalo){
          highestECl2 = 0;
        }
      
        Int_t highestECl3 = -1;
        if (b2ndCalo){
          highestECl3 = 0;
        }
      
        if (!(highestECl > -1 || highestECl2 > -1 || highestECl3 > -1)) continue;
        if(verbosityRH>1){
          if (bRefCalo) cout << "cluster: " << icalo << "\t i:" << highestECl<<"\t" << (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_E << endl;
          if (b1stCalo) cout << "cluster 1: " << _combCalo[icalo] << "\t i:" << highestECl2<<"\t" << (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E << endl;
          if (b2ndCalo) cout << "cluster 2" << _combCalo2[icalo] << "\t i:" << highestECl3<<"\t" << (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E << endl;
        }

        float clustereta                  = 0;
        float clusterenergy               = 0;
        float clusterenergyCalib          = 0;
        float clusterengeryCalib_smeared  = 0;
        float clusterengery_smeared       = 0;
        float tempEcalib                  = 0;
        Int_t currMCID                    = -999999;
        // cluster should have at least 2 towers
        if ( highestECl > -1){
          if((_clusters_calo[ialgo][icalo].at(highestECl)).cluster_NTowers>1 ){
            currMCID                    = (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_trueID;
            clustereta                  = (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_Eta;
            clusterenergy               = (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_E;
            clusterenergyCalib          = (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_E/
                                            getCalibrationValue((_clusters_calo[ialgo][icalo].at(highestECl)).cluster_E, icalo, ialgo);
            clusterengeryCalib_smeared  = clusterenergyCalib*getEnergySmearing(icalo,ialgo);
            clusterengery_smeared       = clusterenergy*getEnergySmearing(icalo,ialgo);
          }
        }
        // find true MC particle for given cluster
        if ( highestECl2 > -1){
          if((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_NTowers>1 ){
            if (bRefCalo && currMCID != -999999){
              if(currMCID==(_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_trueID){
                clusterenergy               += ((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E);
                clusterenergyCalib          += ((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E/
                                                getCalibrationValue((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E, _combCalo[icalo], ialgo));
                tempEcalib                  = ((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E / 
                                                getCalibrationValue((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E, _combCalo[icalo], ialgo));
                clusterengeryCalib_smeared  += tempEcalib*getEnergySmearing(_combCalo[icalo],ialgo);
                clusterengery_smeared       += (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E*getEnergySmearing(icalo,ialgo);
              }
            } else {
              currMCID                      = (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_trueID;
              clustereta                    = (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E;
              clusterenergy                 += ((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E);
              clusterenergyCalib            += ((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E/
                                                getCalibrationValue((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E, _combCalo[icalo], ialgo));
              tempEcalib                    = ((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E / 
                                              getCalibrationValue((_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E, _combCalo[icalo], ialgo));
              clusterengeryCalib_smeared    += tempEcalib*getEnergySmearing(_combCalo[icalo],ialgo);
              clusterengery_smeared         += (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E*getEnergySmearing(icalo,ialgo);
            }
          }
        }
        // if(verbosityRH>1) cout << "\tfound MC:  particle " << currMCID << "\twith E = " << _mcpart_E[currMCID] << " GeV" << endl;
        if ( highestECl3 > -1){
          if((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_NTowers>1 ){
            if ((bRefCalo || b1stCalo) && currMCID != -999999 ){
              if(currMCID==(_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_trueID){
                clusterenergy               += ((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E);
                clusterenergyCalib          += ((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E/
                                              getCalibrationValue((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E, _combCalo2[icalo], ialgo));
                tempEcalib                  = ((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E / 
                                              getCalibrationValue((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E, _combCalo2[icalo], ialgo));
                clusterengeryCalib_smeared  += tempEcalib*getEnergySmearing(_combCalo2[icalo],ialgo);
                clusterengery_smeared       += (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E*getEnergySmearing(icalo,ialgo);
              }
            } else {
              currMCID                      = (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_trueID;
              clustereta                    = (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E;
              clusterenergy                 += ((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E);
              clusterenergyCalib            += ((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E/
                                                getCalibrationValue((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E, _combCalo2[icalo], ialgo));
              tempEcalib                    = ((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E / 
                                                getCalibrationValue((_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E, _combCalo2[icalo], ialgo));
              clusterengeryCalib_smeared    += tempEcalib*getEnergySmearing(_combCalo2[icalo],ialgo);
              clusterengery_smeared         += (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E*getEnergySmearing(icalo,ialgo);
            }
          }
        }
        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[currMCID])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }

        h_RH_ResoComb_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        if (pClass != -1) h_RH_ResoComb_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        h_RH_ResoComb_smear_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterengery_smeared/_mcpart_E[currMCID]);
        if (pClass != -1) h_RH_ResoComb_smear_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterengery_smeared/_mcpart_E[currMCID]);

        Int_t et = 0;
        while ( ( clustereta > partEta[et+1] ) && ( et < nEta )) et++;
        h_RH_ResoComb_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        if (pClass != -1) h_RH_ResoComb_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        
        h_RH_ResoComb_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
        if (pClass != -1) h_RH_ResoComb_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
        
        h_RH_ResoCombCalib_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergyCalib/_mcpart_E[currMCID]);
        if (pClass != -1) h_RH_ResoCombCalib_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergyCalib/_mcpart_E[currMCID]);

        h_RH_ResoCombCalib_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[currMCID]);
        if (pClass != -1) h_RH_ResoCombCalib_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[currMCID]);

        h_RH_ResoCombCalib_smeared_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterengeryCalib_smeared/_mcpart_E[currMCID]);
        if (pClass != -1) h_RH_ResoCombCalib_smeared_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterengeryCalib_smeared/_mcpart_E[currMCID]);
        
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
    if (!caloEnabled[icalo]) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      for(int sp = 0; sp < particleSpRes; sp++){
        if(h_RH_Reso_E[sp][icalo][ialgo])h_RH_Reso_E[sp][icalo][ialgo]->Write();
        if (h_RH_Reso_Ehighest[sp][icalo][ialgo])h_RH_Reso_Ehighest[sp][icalo][ialgo]->Write();
        if(h_RH_Reso_smear_E[sp][icalo][ialgo])h_RH_Reso_smear_E[sp][icalo][ialgo]->Write();
        if(h_RH_Reso_Erec[sp][icalo][ialgo])h_RH_Reso_Erec[sp][icalo][ialgo]->Write();
        if(h_RH_ResoCalib_E[sp][icalo][ialgo])h_RH_ResoCalib_E[sp][icalo][ialgo]->Write();
        if(h_RH_ResoCalib_Erec[sp][icalo][ialgo])h_RH_ResoCalib_Erec[sp][icalo][ialgo]->Write();
        if(h_RH_ResoCalib_smeared_Erec[sp][icalo][ialgo])h_RH_ResoCalib_smeared_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_RH_Reso_E_Eta[sp][icalo][ialgo][eT])h_RH_Reso_E_Eta[sp][icalo][ialgo][eT]->Write();
        }
      }
      
      
      for(int sp = 0; sp < particleSpRes; sp++){
        if(h_RH_ResoComb_E[sp][icalo][ialgo])h_RH_ResoComb_E[sp][icalo][ialgo]->Write();
        if(h_RH_ResoComb_smear_E[sp][icalo][ialgo])h_RH_ResoComb_smear_E[sp][icalo][ialgo]->Write();
        if(h_RH_ResoComb_Erec[sp][icalo][ialgo])h_RH_ResoComb_Erec[sp][icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_E[sp][icalo][ialgo])h_RH_ResoCombCalib_E[sp][icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_Erec[sp][icalo][ialgo])h_RH_ResoCombCalib_Erec[sp][icalo][ialgo]->Write();
        if(h_RH_ResoCombCalib_smeared_Erec[sp][icalo][ialgo])h_RH_ResoCombCalib_smeared_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_RH_ResoComb_E_Eta[sp][icalo][ialgo][eT])h_RH_ResoComb_E_Eta[sp][icalo][ialgo][eT]->Write();
        }
      }
    
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
