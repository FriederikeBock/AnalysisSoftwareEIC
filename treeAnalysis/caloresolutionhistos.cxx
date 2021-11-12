
// ANCHOR debug output verbosity
Int_t verbosityCRH = 0;
const int particleSpRes = 7;
TString addNameRes[particleSpRes]   = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_" ,""};
int pdgCodeRes[particleSpRes]       = {22, 11, 211, 2212, 321, 2112, -1};
bool activeRes[particleSpRes]       = {false, false, false, false, false, false, true};

// ANCHOR create histograms globally
// ========================== energy resol ======================================
TH2F*  h_CRH_EReso_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EReso_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EReso_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];
TH2F*  h_CRH_EReso_Ehighest[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EReso_Ehighest_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];

//========================== ECals and HCal energies together ===================
TH2F*  h_CRH_EResoComb_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoComb_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoComb_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];
TH2F*  h_CRH_EResoCombEMonly_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombEMonly_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombEMonly_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];
TH2F*  h_CRH_EResoCombHighest_E[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombHighest_Erec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_CRH_EResoCombHighest_E_Eta[particleSpRes][_active_calo][_active_algo][nEta+1];

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

struct simpleHighestComb {
  simpleHighestComb(): cluster_MCID(-1), cluster_E(-1000), cluster_Eta(-1000), MC_E(-1000), MC_Eta(-1000), MC_PDG(-10000)  {}
  simpleHighestComb(int mcid, float clE, float clEta, float mcE, float mcEta, int pdg){
    cluster_MCID  = mcid;
    cluster_E     = clE;
    cluster_Eta   = clEta;
    MC_E          = mcE;
    MC_Eta        = mcEta;
    MC_PDG        = pdg;
  }
  int cluster_MCID;
  float cluster_E;
  float cluster_Eta;
  float MC_E;
  float MC_Eta;
  int MC_PDG;
};


// ANCHOR main function to be called in event loop
void caloresolutionhistos(){

  CheckBarrelCaloVersion();
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(verbosityCRH) std::cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << std::endl;
      if(!loadClusterizerInput( ialgo, icalo )) continue;
      
      int nbinsERes = 400;
      int nbinsPhiRes = 200;
      int nbinsEtaRes = 200;
      if(icalo==kFHCAL || icalo==kLFHCAL || icalo==kEEMC || icalo==kBECAL || icalo==kFEMC ){
        nbinsERes = 200;
      } else if(icalo==kEHCAL || icalo==kHCALOUT || icalo==kHCALOUT) {
        nbinsERes = 100;
      }
      
      Int_t caloDir           = GetCaloDirection(icalo);
      Float_t minEtaCurCalo   = partEtaCalo[minEtaBinCalo[caloDir]];
      Float_t maxEtaCurCalo   = partEtaCalo[maxEtaBinCalo[caloDir]];
//       std::cout << "direction: " << caloDir << "\t" << minEtaCurCalo << "\t < \t" << maxEtaCurCalo << std::endl;
      
      for(Int_t sp = 0; sp< particleSpRes; sp++){
        if(!h_CRH_EReso_E[sp][icalo][ialgo]) h_CRH_EReso_E[sp][icalo][ialgo]              = new TH2F(Form("h_CRH_EReso_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EReso_Ehighest[sp][icalo][ialgo]) h_CRH_EReso_Ehighest[sp][icalo][ialgo] = new TH2F(Form("h_CRH_EReso_%sEhighest_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EReso_Erec[sp][icalo][ialgo])h_CRH_EReso_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_CRH_EReso_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
      
        for (Int_t eT = minEtaBinCalo[caloDir]; eT < maxEtaBinCalo[caloDir]; eT++){
//           std::cout << eT << "\t" << partEtaCalo[eT] << "\t < eta < \t" << partEtaCalo[eT+1] << std::endl;
          if(!h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT])h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_CRH_EReso_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
          if(!h_CRH_EReso_Ehighest_Eta[sp][icalo][ialgo][eT])h_CRH_EReso_Ehighest_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_CRH_EReso_%sEhighest_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
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

        // cluster should have at least 2 towers
        if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers<2) continue;
        
        // get MC particle ID
        int mcID = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;
        // loop over MC particles and find particle type
        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[mcID])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
        if (pClass != -1){
          if (!activeRes[pClass]) activeRes[pClass] = true;
        }
        
        bool isHighestForPart = true;
        for(Int_t iclus2=0; iclus2<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus2++){
          if( iclus == iclus2) continue;
          if ( (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_NTowers<2) continue;
          if ( (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID == (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_trueID  ) {
            if ( (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E < (_clusters_calo[ialgo][icalo].at(iclus2)).cluster_E  ){
              isHighestForPart = false;
            }
          }
        }
        if (verbosityCRH > 2) std::cout << icalo << "\t cluster \t" << iclus << "\t eta \t" << (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta << "\t phi \t" << (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi << std::endl;
        
        float mcEta = _mcpart_Eta[mcID];
        float mcPhi = _mcpart_Phi[mcID];
        if (_mcpart_EcalProjs[mcID].size() > 0 && !IsHCALCalorimeter(icalo)){
          int proj = -1;
          for (Int_t i = 0; i < (int)_mcpart_EcalProjs[mcID].size(); i++)
            if (_mcpart_EcalProjs[mcID].at(i).caloid == icalo) proj=i;
          if ( proj != -1){
            mcEta = _mcpart_EcalProjs[mcID].at(proj).eta_Calo;
            mcPhi = _mcpart_EcalProjs[mcID].at(proj).phi_Calo;
            if (verbosityCRH > 2) std::cout << "found projection: " << icalo << "\t"  << mcEta << "\t" << mcPhi << std::endl;
          } else {
            if (verbosityCRH > 2) std::cout << "no correct projection found, using phi & eta at vtx" << std::endl;
          }      
        }
        if (_mcpart_HcalProjs[mcID].size() > 0 && IsHCALCalorimeter(icalo)){
          int proj = -1;
          for (Int_t i = 0; i < (int)_mcpart_HcalProjs[mcID].size(); i++)
            if (_mcpart_HcalProjs[mcID].at(i).caloid == icalo) proj=i;
          if ( proj != -1){
            mcEta = _mcpart_HcalProjs[mcID].at(proj).eta_Calo;
            mcPhi = _mcpart_HcalProjs[mcID].at(proj).phi_Calo;
            if (verbosityCRH > 2) std::cout << "found projection: " << icalo << "\t"  << mcEta << "\t" << mcPhi << std::endl;
          } else {
            if (verbosityCRH > 2) std::cout << "no correct projection found, using phi & eta at vtx" << std::endl;
          }
        }
        
        bool isInNomEta = false;
        if (icalo != kHCALOUT){
          if (mcEta >= nominalEtaRegion[caloDir][0] && mcEta <= nominalEtaRegion[caloDir][1])
            isInNomEta    = true;
        } else {
          if (mcEta >= -1.1 && mcEta <= 1.1)
            isInNomEta    = true;
        }
        // find true MC particle for given cluster
        // if(verbosityCRH>1) std::cout << "\tfound MC:  particle " << mcID << "\twith E = " << _mcpart_E[mcID] << " GeV" << std::endl;
        // ========== Energy resolutions ================================================
        float energyRes = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E/_mcpart_E[mcID];
        // res E vs MC E
        if (isInNomEta){
          h_CRH_EReso_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
          if (pClass != -1) h_CRH_EReso_E[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
          // highest energy cluster res E vs MC E
          if ( isHighestForPart ){
            h_CRH_EReso_Ehighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
            if (pClass != -1) h_CRH_EReso_Ehighest[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],energyRes);
          }
          // res E vs rec E
          h_CRH_EReso_Erec[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,energyRes);
          if (pClass != -1) 
            h_CRH_EReso_Erec[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,energyRes);
        }
        
        // find eta bin
        Int_t et = 0;
        while ( ( _mcpart_Eta[mcID] > partEtaCalo[et+1] ) && ( et < nEta )) et++;
//         while ( ( (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta > partEtaCalo[et+1] ) && ( et < nEta )) et++;
        // res E vs MC E in eta slices
        if(et >= minEtaBinCalo[caloDir] && et < maxEtaBinCalo[caloDir]) {
          h_CRH_EReso_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[mcID],energyRes);
          if (pClass != -1) 
            h_CRH_EReso_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[mcID],energyRes);
          if (isHighestForPart){
            h_CRH_EReso_Ehighest_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[mcID],energyRes);
            if (pClass != -1) 
              h_CRH_EReso_Ehighest_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[mcID],energyRes);            
          }
        }

        //====== Eta resolutions =====================================================//                
        float etadiff = ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-mcEta);
        // res Eta vs MC eta
        h_CRH_EtaReso_Eta[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Eta[mcID], etadiff);
        if (pClass != -1) 
          h_CRH_EtaReso_Eta[pClass][icalo][ialgo]->Fill(_mcpart_Eta[mcID], etadiff);
        // highest energy cluster res Eta vs MC eta
        
        if (isHighestForPart ){
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
        float phidiff = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-mcPhi;
        if (phidiff >= TMath::Pi()) phidiff = phidiff-2*TMath::Pi();
        if (phidiff <= -TMath::Pi()) phidiff = phidiff+2*TMath::Pi();
        
        if (isInNomEta){
          // res Phi vs MC Phi
          h_CRH_PhiReso_Phi[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Phi[mcID], phidiff);
          if (pClass != -1) 
            h_CRH_PhiReso_Phi[pClass][icalo][ialgo]->Fill(_mcpart_Phi[mcID], phidiff);
          // higherst energy cluster res Phi vs MC Phi
          if (isHighestForPart ){
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
  }
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;          // calorimeter is enabled
    if (_combCalo[icalo] == -1 ) continue;      // inner Hcal or Ecal defined
    if (!caloEnabled[_combCalo[icalo]]) continue; // inner HCal or Ecal ran
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if (_combCalo[icalo] == -1 ) continue;
      bool bRefCalo = loadClusterizerInput( ialgo, icalo);
      bool b1stCalo = loadClusterizerInput( ialgo, _combCalo[icalo]);
      Int_t caloDir = GetCaloDirection(icalo);
        
      bool b2ndCalo = false;
      if (_combCalo2[icalo] != -1 ){
        if (caloEnabled[_combCalo2[icalo]]){ 
          b2ndCalo = loadClusterizerInput( ialgo, _combCalo2[icalo]);
        }
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
        if(!h_CRH_EResoComb_Erec[sp][icalo][ialgo])h_CRH_EResoComb_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_CRH_EResoComb_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        for (Int_t eT = minEtaBinCalo[caloDir]; eT < maxEtaBinCalo[caloDir]; eT++){
          if(!h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_CRH_EResoComb_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        }        
        if(!h_CRH_EResoCombEMonly_E[sp][icalo][ialgo]) h_CRH_EResoCombEMonly_E[sp][icalo][ialgo]              = new TH2F(Form("h_CRH_EResoCombEMonly_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EResoCombEMonly_Erec[sp][icalo][ialgo])h_CRH_EResoCombEMonly_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_CRH_EResoCombEMonly_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        for (Int_t eT = minEtaBinCalo[caloDir]; eT < maxEtaBinCalo[caloDir]; eT++){
          if(!h_CRH_EResoCombEMonly_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoCombEMonly_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_CRH_EResoCombEMonly_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        }        
        if(!h_CRH_EResoCombHighest_E[sp][icalo][ialgo]) h_CRH_EResoCombHighest_E[sp][icalo][ialgo]              = new TH2F(Form("h_CRH_EResoCombHighest_%sE_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        if(!h_CRH_EResoCombHighest_Erec[sp][icalo][ialgo])h_CRH_EResoCombHighest_Erec[sp][icalo][ialgo]         = new TH2F(Form("h_CRH_EResoCombHighest_%sErec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        for (Int_t eT = minEtaBinCalo[caloDir]; eT < maxEtaBinCalo[caloDir]; eT++){
          if(!h_CRH_EResoCombHighest_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoCombHighest_E_Eta[sp][icalo][ialgo][eT] 	= new TH2F(Form("h_CRH_EResoCombHighest_%sE_Eta_%s_%s_%d", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data(), eT), "", 500, 0.25, 250.25, nbinsERes, 0, 2);
        }        
      }

      std::vector<simpleHighestComb> highestComb;
      //******************************************************************************
      //Hcal loop
      //******************************************************************************
      for (int iclus = 0; iclus < (int)_clusters_calo[ialgo][icalo].size(); iclus++){
        if(!((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers > 1) ) continue;
        float clustereta                  = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta;
        float clusterenergy               = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E;
        Int_t currMCID                    = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;
        
        
        // find matching ECal cluster
        float dRE   = 10000; 
        int clIDE   = -1;
        int calIDE  = -1;
        float energyE = 0;
        for (int iCl2 = 0; iCl2 < (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals.size()); iCl2++){
          int clID        = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).id;
          int calID       = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).caloid;
          float dPhi      = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).dPhi;
          float dEta      = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedECals).at(iCl2)).dEta;
          float energyCl2 = _clusters_calo[ialgo][calID].at(clID).cluster_E;
          float dRC       = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
          if ( dRC < dRE && energyCl2 > 0.9*energyE && (_clusters_calo[ialgo][calID].at(clID)).cluster_NTowers>1 ){
            if(verbosityCRH >2) std::cout << "\t current\t"<< clIDE << "\t" << calIDE << "\t" << dRE << "\t" << energyE << std::endl;
            if(verbosityCRH >2) std::cout << "\t new\t"<< clID << "\t" << calID << "\t" << dRC << "\t" << energyCl2 << std::endl;
            clIDE   = clID;
            calIDE  = calID;
            energyE = energyCl2;
          }
        }
        if ( clIDE > -1){
          if(currMCID==(_clusters_calo[ialgo][calIDE].at(clIDE)).cluster_trueID){
            clusterenergy               += energyE;
          }
        }

        
        // find matching HCal cluster for barrel kHCALOUT
        float dRH   = 10000; 
        int clIDH   = -1;
        int calIDH  = -1;
        float energyH = 0;
        if (icalo == kHCALOUT){
          for (int iCl3 = 0; iCl3 < (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals.size()); iCl3++){
            int clID        = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl3)).id;
            int calID       = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl3)).caloid;
            float dPhi      = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl3)).dPhi;
            float dEta      = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedHCals).at(iCl3)).dEta;
            float energyCl3 = _clusters_calo[ialgo][calID].at(clID).cluster_E;
            float dRC       = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
            if ( dRC < dRH && energyCl3 > 0.9*energyH && (_clusters_calo[ialgo][calID].at(clID)).cluster_NTowers>1 ){
              if(verbosityCRH >2) std::cout << "\t current\t"<< clIDH << "\t" << calIDH << "\t" << dRH << "\t" << energyH << std::endl;
              if(verbosityCRH >2) std::cout << "\t new\t"<< clID << "\t" << calID << "\t" << dRC << "\t" << energyCl3 << std::endl;
              clIDH   = clID;
              calIDH  = calID;
              energyH = energyCl3;
            }
          }
        }
        
        if (clIDH > -1){
          if(currMCID==(_clusters_calo[ialgo][calIDH].at(clIDH)).cluster_trueID){
            clusterenergy               += energyH;
          }
        }
        simpleHighestComb currentPart(currMCID,clusterenergy,clustereta, _mcpart_E[currMCID], _mcpart_Eta[currMCID], _mcpart_PDG[currMCID]);
        if (!highestComb.empty()){
          bool hasMCpart = false;
          int entry = -1;
          for (int l = 0; l < (int)highestComb.size(); l++){
            if (highestComb.at(l).cluster_MCID == currMCID){
              hasMCpart = true;
              entry     = l;
            }
          }
          if (hasMCpart){
            if (highestComb.at(entry).cluster_E < clusterenergy){
              highestComb.at(entry).cluster_E = clusterenergy;
              highestComb.at(entry).cluster_Eta = clustereta;
            }
          } else {
            highestComb.push_back(currentPart);
          }
        } else {
          highestComb.push_back(currentPart);
        }
      
        if(verbosityCRH>1){
          std::cout << "cluster: " << icalo << "\t i:" << iclus <<"\t" << (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E << std::endl;
          if (b1stCalo) std::cout << "cluster 1: " << _combCalo[icalo] << "\t i:" << clIDE<<"\t" << (_clusters_calo[ialgo][_combCalo[icalo]].at(clIDE)).cluster_E << std::endl;
          if (b2ndCalo) std::cout << "cluster 2" << _combCalo2[icalo] << "\t i:" << clIDH<<"\t" << (_clusters_calo[ialgo][_combCalo2[icalo]].at(clIDH)).cluster_E << std::endl;
        }

        float mcEta                  = _mcpart_Eta[currMCID];
        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[currMCID])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
        
        bool isInNomEta = false;
        if (icalo != kHCALOUT){
          if (mcEta >= nominalEtaRegion[caloDir][0] && mcEta <= nominalEtaRegion[caloDir][1])
            isInNomEta    = true;
        } else {
          if (mcEta >= -1.1 && mcEta <= 1.1)
            isInNomEta    = true;
        }
        
        if (isInNomEta){    
          // combined res E vs MC E
          h_CRH_EResoComb_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
          if (pClass != -1) h_CRH_EResoComb_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
          // combined res E vs rec E
          h_CRH_EResoComb_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
          if (pClass != -1) h_CRH_EResoComb_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
        }
        // find eta bin
        Int_t et = 0;
        while ( ( mcEta > partEtaCalo[et+1] ) && ( et < nEta )) et++;
        // combined res E vs MC E in eta slices
        if(et >= minEtaBinCalo[caloDir] && et < maxEtaBinCalo[caloDir]) {
          h_CRH_EResoComb_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
          if (pClass != -1) h_CRH_EResoComb_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        }
      }

      //******************************************************************************
      // ecal only loop
      //******************************************************************************
      for (int iclus = 0; iclus < (Int_t)_clusters_calo[ialgo][_combCalo[icalo]].size(); iclus++){
        if((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_NTowers<2 ) continue;
          
        bool matchedWHCal = false;
        for (int iCl3 = 0; iCl3 < (int)((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_matchedHCals.size()); iCl3++){
          int clID        = (((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_matchedHCals).at(iCl3)).id;
          int calID       = (((_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_matchedHCals).at(iCl3)).caloid;
          if (_clusters_calo[ialgo][calID].at(clID).cluster_NTowers > 1 )
            matchedWHCal  = true;
        }
        // we only fill here the Ecal cluster which don't have a HCal cluster with Ntowers > 1 matched
        if (matchedWHCal) continue;     
        
        if(verbosityCRH>1){
          std::cout << "cluster Ecal only: " << _combCalo[icalo] << "\t i:" << iclus <<"\t" << (_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_E << std::endl;
        }

        float clustereta                  = (_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_Eta;
        float clusterenergy               = (_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_E;
        Int_t currMCID                    = (_clusters_calo[ialgo][_combCalo[icalo]].at(iclus)).cluster_trueID;

        simpleHighestComb currentPart( currMCID, clusterenergy, clustereta, _mcpart_E[currMCID], _mcpart_Eta[currMCID], _mcpart_PDG[currMCID]);
        if (!highestComb.empty()){
          bool hasMCpart = false;
          int entry = -1;
          for (int l = 0; l < (int)highestComb.size(); l++){
            if (highestComb.at(l).cluster_MCID == currMCID){
              hasMCpart = true;
              entry     = l;
            }
          }
          if (hasMCpart){
            if (highestComb.at(entry).cluster_E < clusterenergy){
              highestComb.at(entry).cluster_E = clusterenergy;
              highestComb.at(entry).cluster_Eta = clustereta;
            }
          } else {
            highestComb.push_back(currentPart);
          }
        } else {
          highestComb.push_back(currentPart);
        }

        
        float mcEta                  = _mcpart_Eta[currMCID];
        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[currMCID])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
        
        bool isInNomEta = false;
        if (icalo != kHCALOUT){
          if (mcEta >= nominalEtaRegion[caloDir][0] && mcEta <= nominalEtaRegion[caloDir][1])
            isInNomEta    = true;
        } else {
          if (mcEta >= -1.1 && mcEta <= 1.1)
            isInNomEta    = true;
        }
        if (isInNomEta){    
          // combined res E vs MC E
          h_CRH_EResoCombEMonly_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
          if (pClass != -1) h_CRH_EResoCombEMonly_E[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
          // combined res E vs rec E
          h_CRH_EResoCombEMonly_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
          if (pClass != -1) h_CRH_EResoCombEMonly_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[currMCID]);
        }
        // find eta bin
        Int_t et = 0;
        while ( ( mcEta > partEtaCalo[et+1] ) && ( et < nEta )) et++;
        // combined res E vs MC E in eta slices
        if(et >= minEtaBinCalo[caloDir] && et < maxEtaBinCalo[caloDir]) {
          h_CRH_EResoCombEMonly_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
          if (pClass != -1) h_CRH_EResoCombEMonly_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[currMCID],clusterenergy/_mcpart_E[currMCID]);
        }
      }
      
      //******************************************************************************
      // highest combined energy for particles
      //******************************************************************************
      for (Int_t iclus =0 ; iclus < (int)highestComb.size(); iclus++){
        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(highestComb.at(iclus).MC_PDG)== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
  
        bool isInNomEta = false;
        if (icalo != kHCALOUT){
          if (highestComb.at(iclus).MC_Eta >= nominalEtaRegion[caloDir][0] && highestComb.at(iclus).MC_Eta <= nominalEtaRegion[caloDir][1])
            isInNomEta    = true;
        } else {
          if (highestComb.at(iclus).MC_Eta >= -1.1 && highestComb.at(iclus).MC_Eta <= 1.1)
            isInNomEta    = true;
        }
        if (isInNomEta){
          // combined res E vs MC E
          h_CRH_EResoCombHighest_E[particleSpRes-1][icalo][ialgo]->Fill(highestComb.at(iclus).MC_E,highestComb.at(iclus).cluster_E/highestComb.at(iclus).MC_E);
          if (pClass != -1) h_CRH_EResoCombHighest_E[pClass][icalo][ialgo]->Fill(highestComb.at(iclus).MC_E,highestComb.at(iclus).cluster_E/highestComb.at(iclus).MC_E);
          // combined res E vs rec E
          h_CRH_EResoCombHighest_Erec[particleSpRes-1][icalo][ialgo]->Fill(highestComb.at(iclus).cluster_E,highestComb.at(iclus).cluster_E/highestComb.at(iclus).MC_E);
          if (pClass != -1) h_CRH_EResoCombHighest_Erec[pClass][icalo][ialgo]->Fill(highestComb.at(iclus).cluster_E,highestComb.at(iclus).cluster_E/highestComb.at(iclus).MC_E);
        }
        // find eta bin
        Int_t et = 0;
        while ( ( highestComb.at(iclus).MC_Eta > partEtaCalo[et+1] ) && ( et < nEta )) et++;
        // combined res E vs MC E in eta slices
        if(et >= minEtaBinCalo[caloDir] && et < maxEtaBinCalo[caloDir]) {
          h_CRH_EResoCombHighest_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(highestComb.at(iclus).MC_E,highestComb.at(iclus).cluster_E/highestComb.at(iclus).MC_E);
          if (pClass != -1) h_CRH_EResoCombHighest_E_Eta[pClass][icalo][ialgo][et]->Fill(highestComb.at(iclus).MC_E,highestComb.at(iclus).cluster_E/highestComb.at(iclus).MC_E);
        }
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
        if (!activeRes[sp]) continue;
        // energy resolutions
        if(h_CRH_EReso_E[sp][icalo][ialgo])h_CRH_EReso_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EReso_Ehighest[sp][icalo][ialgo])h_CRH_EReso_Ehighest[sp][icalo][ialgo]->Write();
        if(h_CRH_EReso_Erec[sp][icalo][ialgo])h_CRH_EReso_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT])h_CRH_EReso_E_Eta[sp][icalo][ialgo][eT]->Write();
          if(h_CRH_EReso_Ehighest_Eta[sp][icalo][ialgo][eT])h_CRH_EReso_Ehighest_Eta[sp][icalo][ialgo][eT]->Write();
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
        if(h_CRH_EResoComb_Erec[sp][icalo][ialgo])h_CRH_EResoComb_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoComb_E_Eta[sp][icalo][ialgo][eT]->Write();
        }
        if(h_CRH_EResoCombEMonly_E[sp][icalo][ialgo])h_CRH_EResoCombEMonly_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCombEMonly_Erec[sp][icalo][ialgo])h_CRH_EResoCombEMonly_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_CRH_EResoCombEMonly_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoCombEMonly_E_Eta[sp][icalo][ialgo][eT]->Write();
        }
        if(h_CRH_EResoCombHighest_E[sp][icalo][ialgo])h_CRH_EResoCombHighest_E[sp][icalo][ialgo]->Write();
        if(h_CRH_EResoCombHighest_Erec[sp][icalo][ialgo])h_CRH_EResoCombHighest_Erec[sp][icalo][ialgo]->Write();
        for (Int_t eT = 0; eT < nEta+1; eT++){
          if(h_CRH_EResoCombHighest_E_Eta[sp][icalo][ialgo][eT])h_CRH_EResoCombHighest_E_Eta[sp][icalo][ialgo][eT]->Write();
        }
      }
    
    }
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
