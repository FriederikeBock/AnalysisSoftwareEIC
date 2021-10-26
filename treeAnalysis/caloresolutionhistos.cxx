
// ANCHOR debug output verbosity
Int_t verbosityCRH = 0;

const int particleSpRes = 7;
TString addNameRes[particleSpRes]   = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_" ,""};
int pdgCodeRes[particleSpRes]       = {22, 11, 211, 2212, 321, 2112, -1};

// ANCHOR create histograms globally
// ========================== energy resol ======================================
TH2F*  h_CRH_EReso_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoSmear_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EReso_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCalib_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCalib_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCalibSmear_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EReso_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];
TH2F*  h_CRH_EReso_Ehighest[particleSpRes][_active_calo][_active_algo];

//========================== ECals and HCal energies together ===================
TH2F*  h_CRH_EResoComb_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombSmear_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoComb_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombCalib_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombCalib_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombCalibSmear_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoComb_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];

// ========================== phi resol ==========================================
TH2F*  h_CRH_EtaReso_Eta[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EtaReso_Etarec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EtaReso_Etahighest[particleSpRes][_active_calo][_active_algo];
TH3F*  h_CRH_EtaReso_E_Eta[particleSpRes][_active_calo][_active_algo];
TH3F*  h_CRH_EtaReso_Ereco_Eta[particleSpRes][_active_calo][_active_algo];

// ========================== phi resol ==========================================
TH2F*  h_CRH_PhiReso_Phi[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_PhiReso_Phirec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_PhiReso_Phihighest[particleSpRes][_active_calo][_active_algo];
TH3F*  h_CRH_PhiReso_E_Phi[particleSpRes][_active_calo][_active_algo];
TH3F*  h_CRH_PhiReso_Ereco_Phi[particleSpRes][_active_calo][_active_algo];


// ANCHOR main function to be called in event loop
void caloresolutionhistos(){

  CheckBarrelCaloVersion();
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(verbosityCRH) std::cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << std::endl;
      if(!loadClusterizerInput( ialgo, icalo )) continue;
      
      std::sort(_clusters_calo[ialgo][icalo].begin(), _clusters_calo[ialgo][icalo].end(), &acompareCl);

      int nbinsERes = 400;
      int nbinsPhiRes = 200;
      int nbinsEtaRes = 200;
      if(icalo==kFHCAL || icalo==kLFHCAL || icalo==kEEMC || icalo==kBECAL || icalo==kFEMC ){
        nbinsERes = 200;
      } else if(icalo==kEHCAL || icalo==kHCALOUT || icalo==kHCALOUT) {
        nbinsERes = 100;
      }
      
      Int_t caloDir           = GetCaloDirection(icalo);
      Float_t minEtaCurCalo   = partEta[minEtaBinFull[caloDir]];
      Float_t maxEtaCurCalo   = partEta[maxEtaBinFull[caloDir]+1];
      
      for(Int_t sp = 0; sp< particleSpRes; sp++){
        if(!h_CRH_EReso_E[sp][icalo][ialgo]) h_CRH_EReso_E[sp][icalo][ialgo]              = new TH2F(Form("h_CRH_EReso_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EReso_Ehighest[sp][icalo][ialgo]) h_CRH_EReso_Ehighest[sp][icalo][ialgo] = new TH2F(Form("h_CRH_EReso_%sEhighest_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EResoSmear_E[sp][icalo][ialgo]) h_CRH_EResoSmear_E[sp][icalo][ialgo]  = new TH2F(Form("h_CRH_EResoSmear_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EReso_Erec[sp][icalo][ialgo])h_CRH_EReso_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_CRH_EReso_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EResoCalib_E[sp][icalo][ialgo])h_CRH_EResoCalib_E[sp][icalo][ialgo]      = new TH2F(Form("h_CRH_EResoCalib_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
      
        if(!h_CRH_EResoCalib_Erec[sp][icalo][ialgo])h_CRH_EResoCalib_Erec[sp][icalo][ialgo] = new TH2F(Form("h_CRH_EResoCalib_%sErec_%s_%s", addNameRes[sp].Data(),  str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EResoCalibSmear_Erec[sp][icalo][ialgo])h_CRH_EResoCalibSmear_Erec[sp][icalo][ialgo] = new TH2F(Form("h_CRH_EResoCalibSmear_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
      
        for (Int_t eT = minEtaBinFull[caloDir]; eT < maxEtaBinFull[caloDir]+1; eT++){
          if(!h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT])h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_CRH_EReso_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        }
        
        if(!h_CRH_EtaReso_Eta[sp][icalo][ialgo]) h_CRH_EtaReso_Eta[sp][icalo][ialgo]                = new TH2F(Form("h_CRH_EtaReso_%sEta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo, nbinsEtaRes, -0.2, 0.2);
        if(!h_CRH_EtaReso_Etahighest[sp][icalo][ialgo]) h_CRH_EtaReso_Etahighest[sp][icalo][ialgo]  = new TH2F(Form("h_CRH_EtaReso_%sEtahighest_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo, nbinsEtaRes, -0.2, 0.2);
        if(!h_CRH_EtaReso_Etarec[sp][icalo][ialgo])h_CRH_EtaReso_Etarec[sp][icalo][ialgo]           = new TH2F(Form("h_CRH_EtaReso_%sEtarec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo, nbinsEtaRes, -0.2, 0.2);
        if(!h_CRH_EtaReso_E_Eta[sp][icalo][ialgo])h_CRH_EtaReso_E_Eta[sp][icalo][ialgo] 	      = new TH3F(Form("h_CRH_EtaReso_%sE_Eta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo,nbinsEtaRes, -0.2, 0.2);
        if(!h_CRH_EtaReso_Ereco_Eta[sp][icalo][ialgo])h_CRH_EtaReso_Ereco_Eta[sp][icalo][ialgo] 	      = new TH3F(Form("h_CRH_EtaReso_%sEreco_Eta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "",500, 0.25, 250.25, (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10), minEtaCurCalo, maxEtaCurCalo, nbinsEtaRes, -0.2, 0.2);

        //====== Phi =====================================================//
        if(!h_CRH_PhiReso_Phi[sp][icalo][ialgo]) h_CRH_PhiReso_Phi[sp][icalo][ialgo]                = new TH2F(Form("h_CRH_PhiReso_%sPhi_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 100, -TMath::Pi(), TMath::Pi(), nbinsPhiRes, -1, 1);
        if(!h_CRH_PhiReso_Phihighest[sp][icalo][ialgo]) h_CRH_PhiReso_Phihighest[sp][icalo][ialgo]  = new TH2F(Form("h_CRH_PhiReso_%sPhihighest_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 100, -TMath::Pi(), TMath::Pi(), nbinsPhiRes, -1, 1);
        if(!h_CRH_PhiReso_Phirec[sp][icalo][ialgo])h_CRH_PhiReso_Phirec[sp][icalo][ialgo]           = new TH2F(Form("h_CRH_PhiReso_%sPhirec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 100,-TMath::Pi(), TMath::Pi(), nbinsPhiRes, -1, 1);
        if(!h_CRH_PhiReso_E_Phi[sp][icalo][ialgo])h_CRH_PhiReso_E_Phi[sp][icalo][ialgo]             = new TH3F(Form("h_CRH_PhiReso_%sE_Phi_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, 100,-TMath::Pi(), TMath::Pi(), nbinsPhiRes, -1, 1);
        if(!h_CRH_PhiReso_Ereco_Phi[sp][icalo][ialgo])h_CRH_PhiReso_Ereco_Phi[sp][icalo][ialgo] 	      = new TH3F(Form("h_CRH_PhiReso_%sEreco_Phi_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, 100,-TMath::Pi(), TMath::Pi(), nbinsPhiRes, -1, 1);
       
        
      }
      
      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus++){
        // if(verbosityCRH>1) std::cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << std::endl;
        // if(verbosityCRH>1) std::cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << std::endl;

        // cluster should have at least 2 towers
        if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers<2) continue;
        
        // get cluster energy and calibrate if necessary
        float clusterE_calo = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E;
        clusterE_calo/=getCalibrationValue((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, icalo, ialgo);
        // get MC particle ID
        int mcID = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;
        // loop over MC particles and find particle type
        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[mcID])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
        
        // find true MC particle for given cluster
        // if(verbosityCRH>1) std::cout << "\tfound MC:  particle " << mcID << "\twith E = " << _mcpart_E[mcID] << " GeV" << std::endl;
        // ========== Energy resolutions ================================================
        float energyRes = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID];
        // res E vs MC E
        h_CRH_EReso_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
        if (pClass != -1) h_CRH_EReso_E[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
        // highest energy cluster res E vs MC E
        if (iclus == 0 && _nMCPart==1 ){
          h_CRH_EReso_Ehighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
          if (pClass != -1) h_CRH_EReso_Ehighest[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
        }
        // res E vs rec E
        h_CRH_EReso_Erec[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,energyRes);
        if (pClass != -1) 
          h_CRH_EReso_Erec[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,energyRes);
        // smeared res E vs MC E 
        float smearedClusE = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E*=getEnergySmearing(icalo,ialgo);
        h_CRH_EResoSmear_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],smearedClusE/_mcpart_E[mcID]);
        if (pClass != -1) 
          h_CRH_EResoSmear_E[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],smearedClusE/_mcpart_E[mcID]);
        //calib res E vs MC E
        h_CRH_EResoCalib_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],clusterE_calo/_mcpart_E[mcID]);
        if (pClass != -1) 
          h_CRH_EResoCalib_E[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],clusterE_calo/_mcpart_E[mcID]);
        //calib res E vs rec E
        h_CRH_EResoCalib_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[mcID]);
        if (pClass != -1) 
          h_CRH_EResoCalib_Erec[pClass][icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[mcID]);
        //smeared calib res E vs rec E
        h_CRH_EResoCalibSmear_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterE_calo,(clusterE_calo*getEnergySmearing(icalo,ialgo))/_mcpart_E[mcID]);
        if (pClass != -1)
          h_CRH_EResoCalibSmear_Erec[pClass][icalo][ialgo]->Fill(clusterE_calo,(clusterE_calo*getEnergySmearing(icalo,ialgo))/_mcpart_E[mcID]);
        // find eta bin
        Int_t et = 0;
        while ( ( (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta > partEta[et+1] ) && ( et < nEta )) et++;
        // res E vs MC E in eta slices
        if(et >= minEtaBinFull[caloDir] && et < maxEtaBinFull[caloDir]+1) {
          h_CRH_EReso_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[mcID],energyRes);
          if (pClass != -1) 
            h_CRH_EReso_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[mcID],energyRes);
        }

        //====== Eta resolutions =====================================================//                
        float etadiff = ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]);
        // res Eta vs MC eta
        h_CRH_EtaReso_Eta[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Eta[mcID], etadiff);
        if (pClass != -1) 
          h_CRH_EtaReso_Eta[pClass][icalo][ialgo]->Fill(_mcpart_Eta[mcID], etadiff);
        // highest energy cluster res Eta vs MC eta
        if (iclus == 0 && _nMCPart==1 ){
          h_CRH_EtaReso_Etahighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Eta[mcID],etadiff);
          if (pClass != -1) 
            h_CRH_EtaReso_Etahighest[pClass][icalo][ialgo]->Fill(_mcpart_Eta[mcID],etadiff);
        }
        // res Eta vs rec eta
        h_CRH_EtaReso_Etarec[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,etadiff);
        if (pClass != -1) 
          h_CRH_EtaReso_Etarec[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta, etadiff);
        // res Eta vs MC E vs MC Eta
        h_CRH_EtaReso_E_Eta[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Eta[mcID], etadiff);
        if (pClass != -1) 
          h_CRH_EtaReso_E_Eta[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Eta[mcID], etadiff);
        // res Eta vs rec E vs rec Eta
        h_CRH_EtaReso_Ereco_Eta[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,
                                                                     (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,
                                                                     etadiff);
        if (pClass != -1) 
          h_CRH_EtaReso_Ereco_Eta[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,
                                                              (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,
                                                              etadiff);
        //====== Phi resolutions =====================================================//
        float phiRecS = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi + TMath::Pi();   // shift rec phi by pi
        float phiMCS  = _mcpart_Phi[mcID] + TMath::Pi();                                      // shift 
        float phidiff = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID];
        if (phidiff >= TMath::Pi()) phidiff = phidiff-2*TMath::Pi();
        if (phidiff <= -TMath::Pi()) phidiff = phidiff+2*TMath::Pi();
        
        // res Phi vs MC Phi
        h_CRH_PhiReso_Phi[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Phi[mcID], phidiff);
        if (pClass != -1) 
          h_CRH_PhiReso_Phi[pClass][icalo][ialgo]->Fill(_mcpart_Phi[mcID], phidiff);
        // higherst energy cluster res Phi vs MC Phi
        if (iclus == 0 && _nMCPart==1 ){
          h_CRH_PhiReso_Phihighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Phi[mcID], phidiff);
          if (pClass != -1) 
            h_CRH_PhiReso_Phihighest[pClass][icalo][ialgo]->Fill(_mcpart_Phi[mcID], phidiff);
        }
        // res Phi vs rec Phi
        h_CRH_PhiReso_Phirec[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi, phidiff);
        if (pClass != -1) 
          h_CRH_PhiReso_Phirec[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi, phidiff);
        // res Phi vs MC E vs MC Phi
        h_CRH_PhiReso_E_Phi[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Phi[mcID], phidiff);
        if (pClass != -1) 
          h_CRH_PhiReso_E_Phi[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Phi[mcID], phidiff);
        // res Phi vs rec E vs rec Phi
        h_CRH_PhiReso_Ereco_Phi[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,
                                                                  (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi,
                                                                  phidiff);
        if (pClass != -1) 
          h_CRH_PhiReso_Ereco_Phi[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,
                                                           ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi),
                                                           phidiff);

          
          
          
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
        Int_t caloDir = GetCaloDirection(icalo);

        std::sort(_clusters_calo[ialgo][_combCalo[icalo]].begin(), _clusters_calo[ialgo][_combCalo[icalo]].end(), &acompareCl);    
          
        bool b2ndCalo = false;
        if (_combCalo2[icalo] != -1 ){
          if (caloEnabled[_combCalo2[icalo]]){ 
            b2ndCalo = loadClusterizerInput( ialgo, _combCalo2[icalo]);
          }
          std::sort(_clusters_calo[ialgo][_combCalo2[icalo]].begin(), _clusters_calo[ialgo][_combCalo2[icalo]].end(), &acompareCl);  
        }
        
        if (!(b1stCalo || b2ndCalo)) continue;
          
        if (verbosityCRH){
          std::cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo ;
          if (b1stCalo)
            std::cout  << " with ECal: "<< _combCalo[icalo] ;
          if (b2ndCalo)
            std::cout << " and 2nd calo: " <<  _combCalo2[icalo] << std::endl;
          else 
            std::cout << std::endl;
        }
        int nbinsERes = 400;
        if(icalo==kFHCAL || icalo==kLFHCAL || icalo==kEEMC || icalo==kBECAL || icalo==kFEMC ) 
          nbinsERes = 200;
        else if(icalo==kEHCAL || icalo==kHCALOUT || icalo==kHCALOUT) 
          nbinsERes = 100;
        
        for(Int_t sp = 0; sp< particleSpRes; sp++){
          if(!h_CRH_EResoComb_E[sp][icalo][ialgo]) h_CRH_EResoComb_E[sp][icalo][ialgo]              = new TH2F(Form("h_CRH_EResoComb_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
          if(!h_CRH_EResoCombSmear_E[sp][icalo][ialgo]) h_CRH_EResoCombSmear_E[sp][icalo][ialgo]  = new TH2F(Form("h_CRH_EResoCombSmear_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
          if(!h_CRH_EResoComb_Erec[sp][icalo][ialgo])h_CRH_EResoComb_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_CRH_EResoComb_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
          if(!h_CRH_EResoCombCalib_E[sp][icalo][ialgo])h_CRH_EResoCombCalib_E[sp][icalo][ialgo]      = new TH2F(Form("h_CRH_EResoCombCalib_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
          if(!h_CRH_EResoCombCalib_Erec[sp][icalo][ialgo])h_CRH_EResoCombCalib_Erec[sp][icalo][ialgo] = new TH2F(Form("h_CRH_EResoCombCalib_%sErec_%s_%s", addNameRes[sp].Data(),  str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
          if(!h_CRH_EResoCombCalibSmear_Erec[sp][icalo][ialgo])h_CRH_EResoCombCalibSmear_Erec[sp][icalo][ialgo] = new TH2F(Form("h_CRH_EResoCombCalibSmear_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
          for (Int_t eT = minEtaBinFull[caloDir]; eT < maxEtaBinFull[caloDir]+1; eT++){
            if(!h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_CRH_EResoComb_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
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
        if(verbosityCRH>1){
          if (bRefCalo) std::cout << "cluster: " << icalo << "\t i:" << highestECl<<"\t" << (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_E << std::endl;
          if (b1stCalo) std::cout << "cluster 1: " << _combCalo[icalo] << "\t i:" << highestECl2<<"\t" << (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E << std::endl;
          if (b2ndCalo) std::cout << "cluster 2" << _combCalo2[icalo] << "\t i:" << highestECl3<<"\t" << (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E << std::endl;
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
        // if(verbosityCRH>1) std::cout << "\tfound MC:  particle " << currMCID << "\twith E = " << _mcpart_E[currMCID] << " GeV" << std::endl;
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
        
        // combined res E vs MC E
        h_CRH_EResoComb_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        if (pClass != -1) h_CRH_EResoComb_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        // combined smeared res E vs MC E
        h_CRH_EResoCombSmear_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterengery_smeared/_mcpart_E[currMCID]);
        if (pClass != -1) h_CRH_EResoCombSmear_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterengery_smeared/_mcpart_E[currMCID]);

        // find eta bin
        Int_t et = 0;
        while ( ( clustereta > partEta[et+1] ) && ( et < nEta )) et++;
        // combined res E vs MC E in eta slices
        if(et >= minEtaBinFull[caloDir] && et < maxEtaBinFull[caloDir]+1) {
          h_CRH_EResoComb_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
          if (pClass != -1) h_CRH_EResoComb_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        }
        // combined res E vs rec E
        h_CRH_EResoComb_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
        if (pClass != -1) h_CRH_EResoComb_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
        // combined calibrated res E vs MC E
        h_CRH_EResoCombCalib_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergyCalib/_mcpart_E[currMCID]);
        if (pClass != -1) h_CRH_EResoCombCalib_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergyCalib/_mcpart_E[currMCID]);
        // combined calibrated res E vs rec E
        h_CRH_EResoCombCalib_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[currMCID]);
        if (pClass != -1) h_CRH_EResoCombCalib_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[currMCID]);
        // combined calibrated smeared res E vs rec E
        h_CRH_EResoCombCalibSmear_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterengeryCalib_smeared/_mcpart_E[currMCID]);
        if (pClass != -1) h_CRH_EResoCombCalibSmear_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterengeryCalib_smeared/_mcpart_E[currMCID]);
        
      }
    }
  }
}
// ANCHOR save function after event loop
void caloresolutionhistossave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_CRH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      for(int sp = 0; sp < particleSpRes; sp++){
        // energy resolutions
        if(h_CRH_EReso_E[sp][icalo][ialgo])h_CRH_EReso_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EReso_Ehighest[sp][icalo][ialgo])h_CRH_EReso_Ehighest[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoSmear_E[sp][icalo][ialgo])h_CRH_EResoSmear_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EReso_Erec[sp][icalo][ialgo])h_CRH_EReso_Erec[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCalib_E[sp][icalo][ialgo])h_CRH_EResoCalib_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCalib_Erec[sp][icalo][ialgo])h_CRH_EResoCalib_Erec[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCalibSmear_Erec[sp][icalo][ialgo])h_CRH_EResoCalibSmear_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT])h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT]->Write();
        }
        // eta resolutions
        if(h_CRH_EtaReso_Eta[sp][icalo][ialgo])h_CRH_EtaReso_Eta[sp][icalo][ialgo]->Write();
        if(h_CRH_EtaReso_Etarec[sp][icalo][ialgo])h_CRH_EtaReso_Etarec[sp][icalo][ialgo]->Write();
        if(h_CRH_EtaReso_Etahighest[sp][icalo][ialgo])h_CRH_EtaReso_Etahighest[sp][icalo][ialgo]->Write();
        if(h_CRH_EtaReso_E_Eta[sp][icalo][ialgo])h_CRH_EtaReso_E_Eta[sp][icalo][ialgo]->Write();
        if(h_CRH_EtaReso_Ereco_Eta[sp][icalo][ialgo])h_CRH_EtaReso_Ereco_Eta[sp][icalo][ialgo]->Write();
        // phi resolutions
        if(h_CRH_PhiReso_Phi[sp][icalo][ialgo])h_CRH_PhiReso_Phi[sp][icalo][ialgo]->Write();
        if(h_CRH_PhiReso_Phirec[sp][icalo][ialgo])h_CRH_PhiReso_Phirec[sp][icalo][ialgo]->Write();
        if(h_CRH_PhiReso_Phihighest[sp][icalo][ialgo])h_CRH_PhiReso_Phihighest[sp][icalo][ialgo]->Write();
        if(h_CRH_PhiReso_E_Phi[sp][icalo][ialgo])h_CRH_PhiReso_E_Phi[sp][icalo][ialgo]->Write();
        if(h_CRH_PhiReso_Ereco_Phi[sp][icalo][ialgo])h_CRH_PhiReso_Ereco_Phi[sp][icalo][ialgo]->Write();
      }
  
      // combined energy resolutions
      for(int sp = 0; sp < particleSpRes; sp++){
        if(h_CRH_EResoComb_E[sp][icalo][ialgo])h_CRH_EResoComb_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCombSmear_E[sp][icalo][ialgo])h_CRH_EResoCombSmear_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoComb_Erec[sp][icalo][ialgo])h_CRH_EResoComb_Erec[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCombCalib_E[sp][icalo][ialgo])h_CRH_EResoCombCalib_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCombCalib_Erec[sp][icalo][ialgo])h_CRH_EResoCombCalib_Erec[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCombCalibSmear_Erec[sp][icalo][ialgo])h_CRH_EResoCombCalibSmear_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT]->Write();
        }
      }
    
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
