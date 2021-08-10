
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

int nclusters_RH = 0;
float* clusters_RH_E            = new float[_maxNclusters];
float* clusters_RH_Eta         = new float[_maxNclusters];
float* clusters_RH_Phi         = new float[_maxNclusters];
float* clusters_RH_M02         = new float[_maxNclusters];
float* clusters_RH_M20         = new float[_maxNclusters];
bool* clusters_RH_isMatched         = new bool[_maxNclusters];
int* clusters_RH_NTower         = new int[_maxNclusters];
int* clusters_RH_trueID       = new int[_maxNclusters];
int* clusters_RH_NtrueID       = new int[_maxNclusters];
float* clusters_RH_X         = new float[_maxNclusters];
float* clusters_RH_Y         = new float[_maxNclusters];
float* clusters_RH_Z         = new float[_maxNclusters];

int nclusters_comb_RH = 0;
float* clusters_comb_RH_E            = new float[_maxNclusters];
float* clusters_comb_RH_Eta         = new float[_maxNclusters];
float* clusters_comb_RH_Phi         = new float[_maxNclusters];
float* clusters_comb_RH_M02         = new float[_maxNclusters];
float* clusters_comb_RH_M20         = new float[_maxNclusters];
bool* clusters_comb_RH_isMatched         = new bool[_maxNclusters];
int* clusters_comb_RH_NTower         = new int[_maxNclusters];
int* clusters_comb_RH_trueID       = new int[_maxNclusters];
int* clusters_comb_RH_NtrueID       = new int[_maxNclusters];
float* clusters_comb_RH_X         = new float[_maxNclusters];
float* clusters_comb_RH_Y         = new float[_maxNclusters];
float* clusters_comb_RH_Z         = new float[_maxNclusters];

int nclusters_comb_RE = 0;
float* clusters_comb_RE_E            = new float[_maxNclusters];
float* clusters_comb_RE_Eta         = new float[_maxNclusters];
float* clusters_comb_RE_Phi         = new float[_maxNclusters];
float* clusters_comb_RE_M02         = new float[_maxNclusters];
float* clusters_comb_RE_M20         = new float[_maxNclusters];
bool* clusters_comb_RE_isMatched         = new bool[_maxNclusters];
int* clusters_comb_RE_NTower         = new int[_maxNclusters];
int* clusters_comb_RE_trueID       = new int[_maxNclusters];
int* clusters_comb_RE_NtrueID       = new int[_maxNclusters];
float* clusters_comb_RE_X         = new float[_maxNclusters];
float* clusters_comb_RE_Y         = new float[_maxNclusters];
float* clusters_comb_RE_Z         = new float[_maxNclusters];

int nclusters_comb_RI = 0;
float* clusters_comb_RI_E            = new float[_maxNclusters];
float* clusters_comb_RI_Eta         = new float[_maxNclusters];
float* clusters_comb_RI_Phi         = new float[_maxNclusters];
float* clusters_comb_RI_M02         = new float[_maxNclusters];
float* clusters_comb_RI_M20         = new float[_maxNclusters];
bool* clusters_comb_RI_isMatched         = new bool[_maxNclusters];
int* clusters_comb_RI_NTower         = new int[_maxNclusters];
int* clusters_comb_RI_trueID       = new int[_maxNclusters];
int* clusters_comb_RI_NtrueID       = new int[_maxNclusters];
float* clusters_comb_RI_X         = new float[_maxNclusters];
float* clusters_comb_RI_Y         = new float[_maxNclusters];
float* clusters_comb_RI_Z         = new float[_maxNclusters];


// ANCHOR main function to be called in event loop
void resolutionhistos(){

  // f_FEMC_calib->SetParameters(paramsFEMCcalib);

  // f_FHCAL_calib->SetParameters(paramsHCALcalib);

  CheckBarrelCaloVersion();
  
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(verbosityRH) cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << endl;
      if(!loadClusterizerInput(
        ialgo,
        icalo,
        nclusters_RH,
        clusters_RH_E,
        clusters_RH_Eta,
        clusters_RH_Phi,
        clusters_RH_M02,
        clusters_RH_M20,
        clusters_RH_isMatched,
        clusters_RH_NTower,
        clusters_RH_trueID,
        clusters_RH_NtrueID,
        clusters_RH_X,
        clusters_RH_Y,
        clusters_RH_Z))continue;
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
      Int_t highestECl = 0;
      Float_t tempclE    = 0; 
      for(Int_t iclus=0; iclus<nclusters_RH; iclus++){
        if(clusters_RH_NTower[iclus]<2) continue;
        if (clusters_RH_E[iclus] > tempclE){
          tempclE = clusters_RH_E[iclus];
          highestECl = iclus;
        }
      }
      for(Int_t iclus=0; iclus<nclusters_RH; iclus++){
        // if(verbosityRH>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _clusters_FHCAL_E[iclus] << " GeV" << endl;
        // if(verbosityRH>1) cout << "\tFHCAL: cluster MC ID " << _clusters_FHCAL_trueID[iclus] << endl;

        // cluster should have at least 2 towers
        if(clusters_RH_NTower[iclus]<2) continue;

        float clusterE_calo = clusters_RH_E[iclus];
        clusterE_calo/=getCalibrationValue(clusters_RH_E[iclus], icalo, ialgo);
        // loop over MC particles

        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[clusters_RH_trueID[iclus]])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
        
        // find true MC particle for given cluster
        // if(verbosityRH>1) cout << "\tfound MC:  particle " << clusters_RH_trueID[iclus] << "\twith E = " << _mcpart_E[clusters_RH_trueID[iclus]] << " GeV" << endl;
        h_RH_Reso_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);
        if (pClass != -1) h_RH_Reso_E[pClass][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);

        if (highestECl == iclus && _nMCPart==1 ){
          h_RH_Reso_Ehighest[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);
          if (pClass != -1) h_RH_Reso_Ehighest[pClass][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);
        
        }
        float smearedClusE = clusters_RH_E[iclus]*=getEnergySmearing(icalo,ialgo);
        h_RH_Reso_smear_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],smearedClusE/_mcpart_E[clusters_RH_trueID[iclus]]);
        if (pClass != -1) h_RH_Reso_smear_E[pClass][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],smearedClusE/_mcpart_E[clusters_RH_trueID[iclus]]);

        Int_t et = 0;
        while ( ( clusters_RH_Eta[iclus] > partEta[et+1] ) && ( et < nEta )) et++;
        h_RH_Reso_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);
        if (pClass != -1) h_RH_Reso_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);
        
        h_RH_Reso_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusters_RH_E[iclus],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);
        if (pClass != -1) h_RH_Reso_Erec[pClass][icalo][ialgo]->Fill(clusters_RH_E[iclus],clusters_RH_E[iclus]/_mcpart_E[clusters_RH_trueID[iclus]]);
        
        h_RH_ResoCalib_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusterE_calo/_mcpart_E[clusters_RH_trueID[iclus]]);
        if (pClass != -1) h_RH_ResoCalib_E[pClass][icalo][ialgo]->Fill(_mcpart_E[clusters_RH_trueID[iclus]],clusterE_calo/_mcpart_E[clusters_RH_trueID[iclus]]);

        h_RH_ResoCalib_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[clusters_RH_trueID[iclus]]);
        if (pClass != -1) h_RH_ResoCalib_Erec[pClass][icalo][ialgo]->Fill(clusterE_calo,clusterE_calo/_mcpart_E[clusters_RH_trueID[iclus]]);

        h_RH_ResoCalib_smeared_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterE_calo,(clusterE_calo*getEnergySmearing(icalo,ialgo))/_mcpart_E[clusters_RH_trueID[iclus]]);
        if (pClass != -1) h_RH_ResoCalib_smeared_Erec[pClass][icalo][ialgo]->Fill(clusterE_calo,(clusterE_calo*getEnergySmearing(icalo,ialgo))/_mcpart_E[clusters_RH_trueID[iclus]]);
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
        if(!loadClusterizerInput(
          ialgo,
          icalo,
          nclusters_comb_RH,
          clusters_comb_RH_E,
          clusters_comb_RH_Eta,
          clusters_comb_RH_Phi,
          clusters_comb_RH_M02,
          clusters_comb_RH_M20,
          clusters_comb_RH_isMatched,
          clusters_comb_RH_NTower,
          clusters_comb_RH_trueID,
          clusters_comb_RH_NtrueID,
          clusters_comb_RH_X,
          clusters_comb_RH_Y,
          clusters_comb_RH_Z))continue;
        if(!loadClusterizerInput(
          ialgo,
          _combCalo[icalo],
          nclusters_comb_RE,
          clusters_comb_RE_E,
          clusters_comb_RE_Eta,
          clusters_comb_RE_Phi,
          clusters_comb_RE_M02,
          clusters_comb_RE_M20,
          clusters_comb_RE_isMatched,
          clusters_comb_RE_NTower,
          clusters_comb_RE_trueID,
          clusters_comb_RE_NtrueID,
          clusters_comb_RE_X,
          clusters_comb_RE_Y,
          clusters_comb_RE_Z))continue;
          
        bool b2ndCalo = false;
        if (_combCalo2[icalo] != -1 ){
          if (caloEnabled[_combCalo2[icalo]]){ 
            b2ndCalo = loadClusterizerInput(
                                                ialgo,
                                                _combCalo2[icalo],
                                                nclusters_comb_RI,
                                                clusters_comb_RI_E,
                                                clusters_comb_RI_Eta,
                                                clusters_comb_RI_Phi,
                                                clusters_comb_RI_M02,
                                                clusters_comb_RI_M20,
                                                clusters_comb_RI_isMatched,
                                                clusters_comb_RI_NTower,
                                                clusters_comb_RI_trueID,
                                                clusters_comb_RI_NtrueID,
                                                clusters_comb_RI_X,
                                                clusters_comb_RI_Y,
                                                clusters_comb_RI_Z
                                            );
          }
        }
        if (verbosityRH){
          cout << "Calo-type: " << icalo << "\t Algorithm: " << ialgo << " with ECal: "<< _combCalo[icalo] ;
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
        Float_t tempclE    = 0; 
        for(Int_t iclus=0; iclus<nclusters_comb_RH; iclus++){
          if(clusters_comb_RH_NTower[iclus]<2) continue;
          if (clusters_comb_RH_E[iclus] > tempclE){
            tempclE = clusters_comb_RH_E[iclus];
            highestECl = iclus;
          }
        }
        Int_t highestECl2 = -1;
        Float_t tempclE2    = 0; 
        for(Int_t iclus=0; iclus<nclusters_comb_RE; iclus++){
          if(clusters_comb_RE_NTower[iclus]<2) continue;
          if (clusters_comb_RE_E[iclus] > tempclE2){
            tempclE2 = clusters_comb_RE_E[iclus];
            highestECl2 = iclus;
          }
        }
        Int_t highestECl3 = -1;
        Float_t tempclE3    = 0;         
        if (b2ndCalo){
          for(Int_t iclus=0; iclus<nclusters_comb_RI; iclus++){
            if(clusters_comb_RI_NTower[iclus]<2) continue;
            if (clusters_comb_RI_E[iclus] > highestECl3){
              tempclE3 = clusters_comb_RI_E[iclus];
              highestECl3 = iclus;
            }
          }
        }
        if (!(highestECl > -1 || highestECl2 > -1)) continue;
        if(verbosityRH>1){
          cout << "cluster: " << icalo << "\t i:" << highestECl<<"\t" << tempclE << endl;
          cout << "cluster 2: " << _combCalo[icalo] << "\t i:" << highestECl2<<"\t" << tempclE2 << endl;
          if (b2ndCalo)cout << "cluster 2" << _combCalo2[icalo] << "\t i:" << highestECl3<<"\t" << tempclE3 << endl;
        }

        float clusterenergy = 0;
        float clusterenergyCalib = 0;
        float clusterengeryCalib_smeared = 0;
        
        // cluster should have at least 2 towers
        if ( highestECl > -1){
          if(clusters_comb_RH_NTower[highestECl]>1 && tempclE > 0){
            clusterenergy = clusters_comb_RH_E[highestECl];
            clusterenergyCalib = clusters_comb_RH_E[highestECl];
            clusterenergyCalib/=getCalibrationValue(clusters_comb_RH_E[highestECl], icalo, ialgo);
            clusterengeryCalib_smeared = clusterenergyCalib*=getEnergySmearing(icalo,ialgo);
          }
        }
        // loop over MC particles

        // find true MC particle for given cluster
        if ( highestECl2 > -1){
          if(clusters_comb_RE_NTower[highestECl2]>1 && tempclE2 > 0){
            if(clusters_comb_RH_trueID[highestECl]==clusters_comb_RE_trueID[highestECl2]){
              clusterenergy += (clusters_comb_RE_E[highestECl2]);
              clusterenergyCalib += (clusters_comb_RE_E[highestECl2]/=getCalibrationValue(clusters_comb_RE_E[highestECl2], _combCalo[icalo], ialgo));
            }
          }
        }
        // if(verbosityRH>1) cout << "\tfound MC:  particle " << clusters_comb_RH_trueID[highestECl] << "\twith E = " << _mcpart_E[clusters_comb_RH_trueID[highestECl]] << " GeV" << endl;
        if (b2ndCalo && highestECl3 > -1){
          if(clusters_comb_RI_NTower[highestECl3]>1  && tempclE3 > 0){
            if(clusters_comb_RH_trueID[highestECl3]==clusters_comb_RI_trueID[highestECl3]){
              clusterenergy += (clusters_comb_RI_E[highestECl3]);
              clusterenergyCalib += (clusters_comb_RI_E[highestECl3]/=getCalibrationValue(clusters_comb_RI_E[highestECl3], _combCalo2[icalo], ialgo));
            }
          }
        }
        int pClass = -1;
        int nTry  = 0;
        while (pClass == -1 && nTry < particleSpRes-1){
          if(abs(_mcpart_PDG[clusters_comb_RH_trueID[highestECl]])== pdgCodeRes[nTry]) pClass = nTry;
          nTry++;  
        }
        h_RH_ResoComb_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],clusterenergy/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        if (pClass != -1) h_RH_ResoComb_E[pClass][icalo][ialgo]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],clusterenergy/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        float smearedClusE = clusters_comb_RH_E[highestECl]*=getEnergySmearing(icalo,ialgo);
        h_RH_ResoComb_smear_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],smearedClusE/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        if (pClass != -1) h_RH_ResoComb_smear_E[pClass][icalo][ialgo]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],smearedClusE/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);

        Int_t et = 0;
        while ( ( clusters_comb_RH_Eta[highestECl] > partEta[et+1] ) && ( et < nEta )) et++;
        h_RH_ResoComb_E_Eta[particleSpRes-1][icalo][ialgo][et]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],clusterenergy/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        if (pClass != -1) h_RH_ResoComb_E_Eta[pClass][icalo][ialgo][et]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],clusterenergy/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        
        h_RH_ResoComb_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        if (pClass != -1) h_RH_ResoComb_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergy/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        
        h_RH_ResoCombCalib_E[particleSpRes-1][icalo][ialgo]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],clusterenergyCalib/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        if (pClass != -1) h_RH_ResoCombCalib_E[pClass][icalo][ialgo]->Fill(_mcpart_E[clusters_comb_RH_trueID[highestECl]],clusterenergyCalib/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);

        h_RH_ResoCombCalib_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        if (pClass != -1) h_RH_ResoCombCalib_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterenergyCalib/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);

        h_RH_ResoCombCalib_smeared_Erec[particleSpRes-1][icalo][ialgo]->Fill(clusterenergy,clusterengeryCalib_smeared/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);
        if (pClass != -1) h_RH_ResoCombCalib_smeared_Erec[pClass][icalo][ialgo]->Fill(clusterenergy,clusterengeryCalib_smeared/_mcpart_E[clusters_comb_RH_trueID[highestECl]]);

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
