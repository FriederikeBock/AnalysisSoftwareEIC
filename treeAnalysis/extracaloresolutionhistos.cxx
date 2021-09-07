
// ANCHOR debug output verbosity
Int_t verbosityERH = 0;


// ANCHOR create histograms globally
//====== Eta =====================================================//
TH2F*  h_ERH_Reso_Eta[particleSpRes][_active_calo][_active_algo];
TH2F*  h_ERH_Reso_Etarec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_ERH_Reso_Etahighest[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_Reso_E_Eta[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_Reso_Ereco_Eta[particleSpRes][_active_calo][_active_algo];

TH2F*  h_ERH_ResoComb_Eta[particleSpRes][_active_calo][_active_algo];
TH2F*  h_ERH_ResoComb_Etarec[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_ResoComb_E_Eta[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_ResoComb_Ereco_Eta[particleSpRes][_active_calo][_active_algo];

//====== Phi =====================================================//

TH2F*  h_ERH_Reso_Phi[particleSpRes][_active_calo][_active_algo];
TH2F*  h_ERH_Reso_Phirec[particleSpRes][_active_calo][_active_algo];
TH2F*  h_ERH_Reso_Phihighest[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_Reso_E_Phi[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_Reso_Ereco_Phi[particleSpRes][_active_calo][_active_algo];

TH2F*  h_ERH_ResoComb_Phi[particleSpRes][_active_calo][_active_algo];
TH2F*  h_ERH_ResoComb_Phirec[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_ResoComb_E_Phi[particleSpRes][_active_calo][_active_algo];
TH3F*  h_ERH_ResoComb_Ereco_Phi[particleSpRes][_active_calo][_active_algo];

// ANCHOR main function to be called in event loop
void extracaloresolutionhistos(){

  CheckBarrelCaloVersion();
  
  for(int icalo=0;icalo<_active_calo;icalo++){
  
    //=====Loop only over enabled calorimeters ===========//
    if (!caloEnabled[icalo]) continue;

    for(int ialgo=0;ialgo<_active_algo;ialgo++){

      if(verbosityERH) cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << endl;

      if(!loadClusterizerInput(ialgo,icalo )) continue;
      
      std::sort(_clusters_calo[ialgo][icalo].begin(), _clusters_calo[ialgo][icalo].end(), &acompareCl);

      
      int nbinsx = 400;
      if(icalo==kFHCAL) nbinsx = 200;
      if(icalo==kEHCAL || icalo==kHCALOUT || icalo==kHCALOUT) nbinsx = 100;
      if(icalo==kLFHCAL) nbinsx = 200;
      if(icalo==kEEMC) nbinsx = 200;
      if(icalo==kBECAL) nbinsx = 200;
      if(icalo==kFEMC) nbinsx = 200;
      
      for(Int_t sp = 0; sp< particleSpRes; sp++){
	//====== Eta =====================================================//
        if(!h_ERH_Reso_Eta[sp][icalo][ialgo]) h_ERH_Reso_Eta[sp][icalo][ialgo]                = new TH2F(Form("h_ERH_Reso_%sEta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, -5, 5, nbinsx, -1, 1);
        if(!h_ERH_Reso_Etahighest[sp][icalo][ialgo]) h_ERH_Reso_Etahighest[sp][icalo][ialgo]  = new TH2F(Form("h_ERH_Reso_%sEtahighest_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, -5, 5, nbinsx, -1, 1);
        if(!h_ERH_Reso_Etarec[sp][icalo][ialgo])h_ERH_Reso_Etarec[sp][icalo][ialgo]           = new TH2F(Form("h_ERH_Reso_%sEtarec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500,-5, 5, nbinsx, -1, 1);
        if(!h_ERH_Reso_E_Eta[sp][icalo][ialgo])h_ERH_Reso_E_Eta[sp][icalo][ialgo] 	      = new TH3F(Form("h_ERH_Reso_%sE_Eta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 150, 500, -5, 5,nbinsx, -1, 1);
        if(!h_ERH_Reso_Ereco_Eta[sp][icalo][ialgo])h_ERH_Reso_Ereco_Eta[sp][icalo][ialgo] 	      = new TH3F(Form("h_ERH_Reso_%sEreco_Eta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 150, 500, -5, 5, nbinsx, -1, 1);
	//====== Phi =====================================================//
	if(!h_ERH_Reso_Phi[sp][icalo][ialgo]) h_ERH_Reso_Phi[sp][icalo][ialgo]                = new TH2F(Form("h_ERH_Reso_%sPhi_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, -5, 5, nbinsx, -1, 1);
        if(!h_ERH_Reso_Phihighest[sp][icalo][ialgo]) h_ERH_Reso_Phihighest[sp][icalo][ialgo]  = new TH2F(Form("h_ERH_Reso_%sPhihighest_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, -5, 5, nbinsx, -1, 1);
        if(!h_ERH_Reso_Phirec[sp][icalo][ialgo])h_ERH_Reso_Phirec[sp][icalo][ialgo]           = new TH2F(Form("h_ERH_Reso_%sPhirec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500,-5, 5, nbinsx, -1, 1);
        if(!h_ERH_Reso_E_Phi[sp][icalo][ialgo])h_ERH_Reso_E_Phi[sp][icalo][ialgo] 	      = new TH3F(Form("h_ERH_Reso_%sE_Phi_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 150, 500,-5, 5, nbinsx, -1, 1);
        if(!h_ERH_Reso_Ereco_Phi[sp][icalo][ialgo])h_ERH_Reso_Ereco_Phi[sp][icalo][ialgo] 	      = new TH3F(Form("h_ERH_Reso_%sEreco_Phi_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0, 150, 500,-5, 5, nbinsx, -1, 1);
     
      } 
      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus++){

        // cluster should have at least 2 towers
        if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers<2) continue;

        float clusterE_calo = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E;
        clusterE_calo/=getCalibrationValue((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E, icalo, ialgo);
        
        int mcID = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;

        // loop over MC particles
        //

        int pClass = -1;
        int nTry  = 0;
        
	while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[mcID])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
       
	cout << _mcpart_Eta[mcID] << "  " <<  ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]) << endl;
	//====== Eta =====================================================//        
	h_ERH_Reso_Eta[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Eta[mcID],((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]));
	if (pClass != -1) h_ERH_Reso_Eta[pClass][icalo][ialgo]->Fill(_mcpart_Eta[mcID],((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]));
	//cout << ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]) << endl;
        if (iclus == 0 && _nMCPart==1 ){
	  h_ERH_Reso_Etahighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Eta[mcID],((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]));
	  if (pClass != -1) h_ERH_Reso_Etahighest[pClass][icalo][ialgo]->Fill(_mcpart_Eta[mcID],((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]));
        }
	h_ERH_Reso_Etarec[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]));
        if (pClass != -1) h_ERH_Reso_Etarec[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]));
	
	//====== Phi =====================================================//
	h_ERH_Reso_Phi[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Phi[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);
        if (pClass != -1) h_ERH_Reso_Phi[pClass][icalo][ialgo]->Fill(_mcpart_Phi[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);
        if (iclus == 0 && _nMCPart==1 ){
          h_ERH_Reso_Phihighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Phi[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);
          if (pClass != -1) h_ERH_Reso_Phihighest[pClass][icalo][ialgo]->Fill(_mcpart_Phi[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);
        }
        h_ERH_Reso_Phirec[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);
        if (pClass != -1) h_ERH_Reso_Phirec[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);

	//====== Eta =====================================================//  
	h_ERH_Reso_E_Eta[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Eta[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]);
        if (pClass != -1) h_ERH_Reso_E_Eta[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Eta[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]);

        h_ERH_Reso_Ereco_Eta[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]);
        if (pClass != -1) h_ERH_Reso_Ereco_Eta[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta-_mcpart_Eta[mcID]);
	
	//====== Phi =====================================================//
	h_ERH_Reso_E_Phi[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Phi[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);
        if (pClass != -1) h_ERH_Reso_E_Phi[pClass][icalo][ialgo]->Fill(_mcpart_E[mcID],_mcpart_Phi[mcID],(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);

         h_ERH_Reso_Ereco_Phi[particleSpRes-1][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi,(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);
         if (pClass != -1) h_ERH_Reso_Ereco_Phi[pClass][icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_E,((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi),(_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi-_mcpart_Phi[mcID]);

      
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
          
        if (verbosityERH){
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
     
	  if(!h_ERH_ResoComb_Eta[sp][icalo][ialgo]) h_ERH_ResoComb_Eta[sp][icalo][ialgo]              = new TH2F(Form("h_ERH_ResoComb_%sEta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, -5, 5, nbinsx, -2, 2);
          if(!h_ERH_ResoComb_Etarec[sp][icalo][ialgo])h_ERH_ResoComb_Etarec[sp][icalo][ialgo]         = new TH2F(Form("h_ERH_ResoComb_%sEtarec_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, -5, 5, nbinsx, -2, 2);
	  if(!h_ERH_ResoComb_E_Eta[sp][icalo][ialgo])h_ERH_ResoComb_E_Eta[sp][icalo][ialgo]   = new TH3F(Form("h_ERH_ResoComb_%sE_Eta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0., 150, 500, -5, 5, nbinsx, -2, 2);
          if(!h_ERH_ResoComb_Ereco_Eta[sp][icalo][ialgo])h_ERH_ResoComb_Ereco_Eta[sp][icalo][ialgo]   = new TH3F(Form("h_ERH_ResoComb_%sEreco_Eta_%s_%s", addNameRes[sp].Data(), str_calorimeter[icalo].Data(), str_clusterizer[ialgo].Data()), "", 500, 0., 150, 500, -5, 5, nbinsx, -2, 2);       

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
        if(verbosityERH>1){
          if (bRefCalo) cout << "cluster: " << icalo << "\t i:" << highestECl<<"\t" << (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_E << endl;
          if (b1stCalo) cout << "cluster 1: " << _combCalo[icalo] << "\t i:" << highestECl2<<"\t" << (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_E << endl;
          if (b2ndCalo) cout << "cluster 2" << _combCalo2[icalo] << "\t i:" << highestECl3<<"\t" << (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_E << endl;
        }

        float clustereta                  = 0;
        float clusterphi                  = 0;
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
            clusterphi                  = (_clusters_calo[ialgo][icalo].at(highestECl)).cluster_Phi;
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
              clustereta                    = (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_Eta;
              clusterphi                    = (_clusters_calo[ialgo][_combCalo[icalo]].at(highestECl2)).cluster_Phi;
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
        // if(verbosityERH>1) cout << "\tfound MC:  particle " << currMCID << "\twith E = " << _mcpart_E[currMCID] << " GeV" << endl;
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
              clustereta                    = (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_Eta;
              clusterphi                    = (_clusters_calo[ialgo][_combCalo2[icalo]].at(highestECl3)).cluster_Phi;
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
	
       	//====== Eta =====================================================// 
	h_ERH_ResoComb_Eta[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Eta[currMCID],(clustereta-_mcpart_Eta[currMCID]));
	if (pClass != -1) h_ERH_ResoComb_Eta[pClass][icalo][ialgo]->Fill(_mcpart_Eta[currMCID],(clustereta-_mcpart_Eta[currMCID]));
	h_ERH_ResoComb_Etarec[particleSpRes-1][icalo][ialgo]->Fill(clustereta,(clustereta-_mcpart_Eta[currMCID]));
	if (pClass != -1) h_ERH_ResoComb_Etarec[pClass][icalo][ialgo]->Fill(clustereta,(clustereta-_mcpart_Eta[currMCID]));
	//====== Phi =====================================================// 
	h_ERH_ResoComb_Phi[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_Phi[currMCID],(clusterphi-_mcpart_Phi[currMCID]));
        if (pClass != -1) h_ERH_ResoComb_Phi[pClass][icalo][ialgo]->Fill(_mcpart_Phi[currMCID],(clusterphi-_mcpart_Phi[currMCID]));
        h_ERH_ResoComb_Phirec[particleSpRes-1][icalo][ialgo]->Fill(clustereta,(clusterphi -_mcpart_Phi[currMCID]));
        if (pClass != -1) h_ERH_ResoComb_Phirec[pClass][icalo][ialgo]->Fill(clustereta,(clusterphi - _mcpart_Phi[currMCID]));

	//====== Eta =====================================================// 	
	 h_ERH_ResoComb_E_Eta[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],_mcpart_Eta[currMCID],clustereta-_mcpart_Eta[currMCID]);
        if (pClass != -1) h_ERH_ResoComb_E_Eta[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],_mcpart_Eta[currMCID],clustereta-_mcpart_Eta[currMCID]);
        h_ERH_ResoComb_Ereco_Eta[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,_mcpart_Eta[currMCID],clustereta-_mcpart_Eta[currMCID]);
        if (pClass != -1) h_ERH_ResoComb_Ereco_Eta[pClass][icalo][ialgo]->Fill(clustereta,_mcpart_Eta[currMCID],clustereta-_mcpart_Eta[currMCID]);
	//====== Phi =====================================================// 
	h_ERH_ResoComb_E_Phi[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[currMCID],_mcpart_Phi[currMCID],clusterphi-_mcpart_Phi[currMCID]);
        if (pClass != -1) h_ERH_ResoComb_E_Phi[pClass][icalo][ialgo]->Fill(_mcpart_E[currMCID],_mcpart_Phi[currMCID],clusterphi-_mcpart_Phi[currMCID]);
        h_ERH_ResoComb_Ereco_Phi[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterphi,clusterphi-_mcpart_Phi[currMCID]);
        if (pClass != -1) h_ERH_ResoComb_Ereco_Phi[pClass][icalo][ialgo]->Fill(clustereta,clusterphi,clusterphi-_mcpart_Phi[currMCID]);	

     }
    }
  }
}
// ANCHOR save function after event loop

void extracaloresolutionhistosSave(){
  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_ERH.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      //====== Eta =====================================================// 
      for(int sp = 0; sp < particleSpRes; sp++){
        if(h_ERH_Reso_Eta[sp][icalo][ialgo])h_ERH_Reso_Eta[sp][icalo][ialgo]->Write();
        if(h_ERH_Reso_Etarec[sp][icalo][ialgo])h_ERH_Reso_Etarec[sp][icalo][ialgo]->Write();
        if (h_ERH_Reso_Etahighest[sp][icalo][ialgo])h_ERH_Reso_Etahighest[sp][icalo][ialgo]->Write();
	if(h_ERH_Reso_E_Eta[sp][icalo][ialgo])h_ERH_Reso_E_Eta[sp][icalo][ialgo]->Write();
        if(h_ERH_Reso_Ereco_Eta[sp][icalo][ialgo])h_ERH_Reso_Ereco_Eta[sp][icalo][ialgo]->Write();
      } 
     
      //====== Phi =====================================================// 
      for(int sp = 0; sp < particleSpRes; sp++){
        if(h_ERH_Reso_Phi[sp][icalo][ialgo])h_ERH_Reso_Phi[sp][icalo][ialgo]->Write();
        if(h_ERH_Reso_Phirec[sp][icalo][ialgo])h_ERH_Reso_Phirec[sp][icalo][ialgo]->Write();
        if (h_ERH_Reso_Phihighest[sp][icalo][ialgo])h_ERH_Reso_Phihighest[sp][icalo][ialgo]->Write();
	if(h_ERH_Reso_E_Phi[sp][icalo][ialgo])h_ERH_Reso_E_Phi[sp][icalo][ialgo]->Write();
        if(h_ERH_Reso_Ereco_Phi[sp][icalo][ialgo])h_ERH_Reso_Ereco_Phi[sp][icalo][ialgo]->Write();
      }
    }	
  }

  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
