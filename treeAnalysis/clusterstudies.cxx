
// ANCHOR debug output verbosity
int int_CS_verbosity = 0;
bool createdHistosCLS = false;

const int nParticlesPDG_CS = 7;
int int_CS_mcparticles_PDG[nParticlesPDG_CS] = {11 /*e*/,22 /*gamma*/,211  /*pi*/,  2212/*p*/,  2112/*n*/,  321/*K*/, 13/*mu*/};
TString str_CS_mcparticles[nParticlesPDG_CS] = {"electron", "photon", "cpion", "proton", "neutron", "ckaon", "muon"};

// ANCHOR create histograms globally
TH1F*  h_CS_clusters_E[_active_calo][_active_algo] = {{NULL}};
TH2F*  h_CS_clusters_M02_E[_active_calo][_active_algo] = {{NULL}};
TH2F*  h_CS_clusters_M02_part_E_truth[_active_calo][_active_algo][nParticlesPDG_CS] = {{{NULL}}};
TH2F*  h_CS_clusters_M20_E[_active_calo][_active_algo] = {{NULL}};
TH2F*  h_CS_clusters_M20_part_E_truth[_active_calo][_active_algo][nParticlesPDG_CS] = {{{NULL}}};
TH2F*  h_CS_clusters_NTower_E[_active_calo][_active_algo] = {{NULL}};
TH1D*  h_CS_clusters_NTowerMean_E[_active_calo][_active_algo] = {{NULL}};

TH2F*  h_clusterizer_clsspec_matched_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecMC_matched_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_clusterizer_clsspec_PDG[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_clusterizer_clsspec_matched_PDG[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]

TH2F*  h_clusterizer_1tower_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_1towerMC_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspec_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecMC_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspec_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecMC_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspec_matched_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecMC_matched_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enu
TH2F*  h_clusterizer_clsspec_matched_particle_E_EoverP[_active_calo][_active_algo][nParticlesPDG_CS][nEta]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecMC_matched_particle_E_EoverP[_active_calo][_active_algo][nParticlesPDG_CS][nEta]; // [calorimeter_enum][algorithm_enu

TH2F*  h_clusterizer_clsspecMC_particle_eta_phi[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspec_particle_eta_phi[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSE_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSEMC_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSE_matched_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSEMC_matched_E_eta[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSE_E_NCl[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSEMC_E_NCl[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH1D*  h_clusterizer_NClMean_E[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH1D*  h_clusterizer_NClMean_MCE[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]

TH1F*  h_clusterizer_NClPerMCLabel[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_clusterizer_ClSameMCLabel_dR[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]

TH2F*  h_clusterizer_clsspecSE_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSEMC_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSE_matched_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSEMC_matched_particle_E_eta[_active_calo][_active_algo][nParticlesPDG_CS]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSE_matched_particle_E_EoverP[_active_calo][_active_algo][nParticlesPDG_CS][nEta]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[_active_calo][_active_algo][nParticlesPDG_CS][nEta]; // [calorimeter_enum][algorithm_enum]

//*************************************************************************************************
//*************************************************************************************************
// ANCHOR main function to be called in event loop
//*************************************************************************************************
//*************************************************************************************************
void clusterstudies(){
  //*************************************************************************************************
  // initialize histos if not yet done
  //*************************************************************************************************        
  if (!createdHistosCLS){
    for(int icalo=0;icalo<_active_calo;icalo++){
      if (!caloEnabled[icalo]) continue;
      for(int ialgo=0;ialgo<_active_algo;ialgo++){
        Int_t caloDir           = GetCaloDirection(icalo);
        Float_t minEtaCurCalo   = partEtaCalo[minEtaBinCalo[caloDir]];
        Float_t maxEtaCurCalo   = partEtaCalo[maxEtaBinCalo[caloDir]];
        Int_t nBinsEta          = (Int_t)(TMath::Abs(minEtaCurCalo-maxEtaCurCalo)*10);
        
        int maxNTowers  = 100;
        if (icalo == kLFHCAL ) maxNTowers = 100*5;
        float maxM02    = 4;
        int maxNM02     = 400;
        if (IsHCALCalorimeter(icalo)){
          maxM02      = 5;
          maxNM02     = 250;
        }
      
        h_CS_clusters_E[icalo][ialgo]                 = new TH1F(Form("h_CS_clusters_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP);
        h_CS_clusters_M02_E[icalo][ialgo]             = new TH2F(Form("h_CS_clusters_M02_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "" , maxNM02, 0, maxM02, nBinsP, binningP);
        h_CS_clusters_M20_E[icalo][ialgo]             = new TH2F(Form("h_CS_clusters_M20_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "" , maxNM02, 0, maxM02, nBinsP, binningP);
        h_CS_clusters_NTower_E[icalo][ialgo]          = new TH2F(Form("h_CS_clusters_NTower_E_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "" , nBinsP, binningP, maxNTowers, -0.5, maxNTowers-0.5);

        h_clusterizer_1tower_E_eta[icalo][ialgo]      = new TH2F(Form("h_clusterizer_1tower_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
        h_clusterizer_1towerMC_E_eta[icalo][ialgo]    = new TH2F(Form("h_clusterizer_1towerMC_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);

        h_clusterizer_clsspec_E_eta[icalo][ialgo]     = new TH2F(Form("h_clusterizer_clsspec_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
        h_clusterizer_clsspecMC_E_eta[icalo][ialgo]   = new TH2F(Form("h_clusterizer_clsspecMC_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);

        h_clusterizer_clsspecSE_E_eta[icalo][ialgo]   = new TH2F(Form("h_clusterizer_clsspecSE_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
        h_clusterizer_clsspecSEMC_E_eta[icalo][ialgo] = new TH2F(Form("h_clusterizer_clsspecSEMC_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
        
        if (tracksEnabled){
          h_clusterizer_clsspec_matched_E_eta[icalo][ialgo]     = new TH2F(Form("h_clusterizer_clsspec_matched_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          h_clusterizer_clsspecMC_matched_E_eta[icalo][ialgo]   = new TH2F(Form("h_clusterizer_clsspecMC_matched_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          h_clusterizer_clsspecSE_matched_E_eta[icalo][ialgo]   = new TH2F(Form("h_clusterizer_clsspecSE_matched_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          h_clusterizer_clsspecSEMC_matched_E_eta[icalo][ialgo] = new TH2F(Form("h_clusterizer_clsspecSEMC_matched_E_eta_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          h_clusterizer_clsspec_matched_PDG[icalo][ialgo]       = new TH1F(Form("h_clusterizer_clsspec_matched_PDG_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 2500,0,2500);
        }
        
        h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]     = new TH2F(Form("h_clusterizer_clsspecSE_E_NCl_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, 80, 0.5, 80.5);
        h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]   = new TH2F(Form("h_clusterizer_clsspecSEMC_E_NCl_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", nBinsP, binningP, 80, 0.5, 80.5);

        h_clusterizer_clsspec_PDG[icalo][ialgo]         = new TH1F(Form("h_clusterizer_clsspec_PDG_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 2500,0,2500);

        h_clusterizer_NClPerMCLabel[icalo][ialgo]       = new TH1F(Form("h_clusterizer_NClPerMCLabel_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 10,-0.5,9.5);
        h_clusterizer_ClSameMCLabel_dR[icalo][ialgo]    = new TH1F(Form("h_clusterizer_ClSameMCLabel_dR_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "",600,0,600);

        for (int iPDG = 0; iPDG < nParticlesPDG_CS; iPDG++){
          h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]            = new TH2F(Form("h_CS_clusters_M02_part_E_truth_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "" , maxNM02, 0, maxM02, nBinsP, binningP);
          h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]            = new TH2F(Form("h_CS_clusters_M20_part_E_truth_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "" , maxNM02, 0, maxM02, nBinsP, binningP);
          h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][iPDG]      = new TH2F(Form("h_clusterizer_clsspec_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][iPDG]    = new TH2F(Form("h_clusterizer_clsspecMC_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          h_clusterizer_clsspecMC_particle_eta_phi[icalo][ialgo][iPDG]  = new TH2F(Form("h_clusterizer_clsspecMC_particle_eta_phi_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsEta*4, minEtaCurCalo,maxEtaCurCalo, 400, -TMath::Pi(), TMath::Pi());
          h_clusterizer_clsspec_particle_eta_phi[icalo][ialgo][iPDG]  = new TH2F(Form("h_clusterizer_clsspec_particle_eta_phi_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsEta*4, minEtaCurCalo,maxEtaCurCalo, 400, -TMath::Pi(), TMath::Pi());
          h_clusterizer_clsspecSE_particle_E_eta[icalo][ialgo][iPDG]    = new TH2F(Form("h_clusterizer_clsspecSE_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          h_clusterizer_clsspecSEMC_particle_E_eta[icalo][ialgo][iPDG]  = new TH2F(Form("h_clusterizer_clsspecSEMC_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
          if (tracksEnabled){
            h_clusterizer_clsspec_matched_particle_E_eta[icalo][ialgo][iPDG]      = new TH2F(Form("h_clusterizer_clsspec_matched_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
            h_clusterizer_clsspecMC_matched_particle_E_eta[icalo][ialgo][iPDG]    = new TH2F(Form("h_clusterizer_clsspecMC_matched_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
            h_clusterizer_clsspecSE_matched_particle_E_eta[icalo][ialgo][iPDG]    = new TH2F(Form("h_clusterizer_clsspecSE_matched_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
            h_clusterizer_clsspecSEMC_matched_particle_E_eta[icalo][ialgo][iPDG]  = new TH2F(Form("h_clusterizer_clsspecSEMC_matched_particle_E_eta_%s_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, nBinsEta, minEtaCurCalo,maxEtaCurCalo);
            
            for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
              Double_t etaMin = partEtaCalo[0];
              Double_t etaMax = partEtaCalo[nEta];
              if (et < nEta){
                etaMin = partEtaCalo[et];
                etaMax = partEtaCalo[et+1];
              }
              h_clusterizer_clsspec_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]      = new TH2F(Form("h_clusterizer_clsspec_matched_particle_Eta_%1.1f_%1.1f_E_EoverP_%s_%s_%s",etaMin,etaMax,str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, 200, 0, 2 );
              h_clusterizer_clsspecMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]    = new TH2F(Form("h_clusterizer_clsspecMC_matched_particle_Eta_%1.1f_%1.1f_E_EoverP_%s_%s_%s",etaMin,etaMax,str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, 200, 0, 2 );
              h_clusterizer_clsspecSE_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]    = new TH2F(Form("h_clusterizer_clsspecSE_matched_particle_Eta_%1.1f_%1.1f_E_EoverP_%s_%s_%s",etaMin,etaMax,str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, 200, 0, 2 );
              h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]  = new TH2F(Form("h_clusterizer_clsspecSEMC_matched_particle_Eta_%1.1f_%1.1f_E_EoverP_%s_%s_%s",etaMin,etaMax,str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data(),str_CS_mcparticles[iPDG].Data()), "", nBinsP, binningP, 200, 0, 2 );
            }
          }
        }
      }
    }
    createdHistosCLS = true;
  }
     
  //*************************************************************************************************
  // fill cluster QA studies for all clusterizer and algorithms
  //*************************************************************************************************         
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){

      if(!loadClusterizerInput(
        ialgo,
        icalo
      )) continue;
      
      std::vector<occuranceStrct> uniqueIDs;
      std::vector<float> usedIDs;
      for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[ialgo][icalo].size(); iclus++){
        Bool_t particleFilled   = kFALSE;
        Double_t highestEnergy  = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E;
        Double_t currEnergy     = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_E;
        Bool_t isHighest        = kFALSE;
        Int_t nClPerParticle    = 1;
        int mcIDCurr            = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_trueID;
        Double_t currEnergyMC   = _mcpart_E[mcIDCurr];
        for (Int_t bclus = iclus+1; bclus < (Int_t)_clusters_calo[ialgo][icalo].size(); bclus++ ){
            if (mcIDCurr != (_clusters_calo[ialgo][icalo].at(bclus)).cluster_trueID) continue;
            if ((_clusters_calo[ialgo][icalo].at(bclus)).cluster_E > currEnergy)
              highestEnergy = (_clusters_calo[ialgo][icalo].at(bclus)).cluster_E;
            nClPerParticle++;
        }
      
        if (currEnergy == highestEnergy)
          isHighest = kTRUE;
        for (int id = 0; id < (int)uniqueIDs.size() && !particleFilled ; id++){
          if ( uniqueIDs.at(id).particle_ID == (int)(mcIDCurr) )
            particleFilled = kTRUE;
        }
        
        usedIDs.push_back(iclus);
        int icluswithsamelabel = 1;
        int icluswithsamelabelM02 = 0;
        if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers>1) icluswithsamelabelM02++;
        for(Int_t icluslbl=0; icluslbl<(Int_t)_clusters_calo[ialgo][icalo].size(); icluslbl++){
          if((int)(mcIDCurr)==(int)((_clusters_calo[ialgo][icalo].at(icluslbl)).cluster_trueID)){
            if ( !(std::find(usedIDs.begin(), usedIDs.end(), icluslbl) != usedIDs.end()) ){
              usedIDs.push_back(icluslbl);
              icluswithsamelabel++;
              h_clusterizer_ClSameMCLabel_dR[icalo][ialgo]->Fill(TMath::Sqrt(TMath::Power((_clusters_calo[ialgo][icalo].at(iclus)).cluster_X-(_clusters_calo[ialgo][icalo].at(icluslbl)).cluster_X,2)+TMath::Power((_clusters_calo[ialgo][icalo].at(iclus)).cluster_Y-(_clusters_calo[ialgo][icalo].at(icluslbl)).cluster_Y,2)));
              if((_clusters_calo[ialgo][icalo].at(icluslbl)).cluster_NTowers>1)icluswithsamelabelM02++;
            }
          }
        }
        
        Double_t etaCurr    = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Eta;
        Double_t phiCurr    = (_clusters_calo[ialgo][icalo].at(iclus)).cluster_Phi;
        Double_t etaCurrMC  = _mcpart_Eta[mcIDCurr];
        
        // determine eta bin
        Int_t et = 0;
        while (partEtaCalo[et+1] < etaCurrMC && et < nEta) et++;

        Double_t trackP = 0;
        if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
          int trackID = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(0)).id;
          TVector3 recpartvec(_track_px[trackID],_track_py[trackID],_track_pz[trackID]);
          int caloIDMatchedTrack  = _track_matchECal[trackID].caloid; 
          int clIDMatchedTrack    = _track_matchECal[trackID].id;
          if (IsHCALCalorimeter(icalo)){
            caloIDMatchedTrack  = _track_matchHCal[trackID].caloid; 
            clIDMatchedTrack    = _track_matchHCal[trackID].id;
          }
          
          trackP = recpartvec.Mag();
          if ( (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).size() > 1){
            for (Int_t k = 1; k < (int)((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).size(); k++){
              int trackIDAlt = (((_clusters_calo[ialgo][icalo].at(iclus)).cluster_matchedTracks).at(k)).id;
              TVector3 recpartvec2(_track_px[trackIDAlt],_track_py[trackIDAlt],_track_pz[trackIDAlt]);
              if (int_CS_verbosity){
                if (k == 1){
                  std::cout << "calo: "<< str_calorimeter[icalo].Data()<< ", found more than one track for cluster " << iclus << "\t" << currEnergy << "\t eta \t"<< etaCurr<< "\t phi \t"<< phiCurr<< std::endl;
                  std::cout << "\t\t\t\t 0\t track P\t" << trackP << "\t eta \t"<< recpartvec.Eta()<< "\t phi \t"<< recpartvec.Phi() << "\t source \t" <<  _track_source[trackID] << "\t calo-enum\t "<<caloIDMatchedTrack << "\t cl ID\t" << clIDMatchedTrack << std::endl; 
                } 
                std::cout << "\t\t\t\t " << k << "\t track P\t" << recpartvec2.Mag() << "\t eta \t"<< recpartvec2.Eta()<< "\t phi \t"<< recpartvec2.Phi()<< "\t source \t" <<  _track_source[trackIDAlt]<< std::endl;
              }
              if (recpartvec2.Mag() > trackP)
                trackP = recpartvec2.Mag();
            }
          }
        }
        
        h_clusterizer_NClPerMCLabel[icalo][ialgo]->Fill(icluswithsamelabel);
        h_CS_clusters_NTower_E[icalo][ialgo]->Fill(currEnergy, (_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers);
        h_CS_clusters_E[icalo][ialgo]->Fill(currEnergy);
        h_CS_clusters_M02_E[icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M02,currEnergy);
        h_CS_clusters_M20_E[icalo][ialgo]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M20,currEnergy);

        for (int iPDG = 0; iPDG < nParticlesPDG_CS; iPDG++){
          if(abs(_mcpart_PDG[mcIDCurr]) == abs(int_CS_mcparticles_PDG[iPDG])){
            h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M02,currEnergy);
            h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]->Fill((_clusters_calo[ialgo][icalo].at(iclus)).cluster_M20,currEnergy);
          }
        }
        if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_NTowers < 2){
          h_clusterizer_1tower_E_eta[icalo][ialgo]->Fill(currEnergy, etaCurr);
          h_clusterizer_1towerMC_E_eta[icalo][ialgo]->Fill(currEnergyMC, etaCurrMC);
          if (!particleFilled && isHighest){
            h_clusterizer_clsspecSE_E_eta[icalo][ialgo]->Fill(currEnergy, etaCurr);
            h_clusterizer_clsspecSEMC_E_eta[icalo][ialgo]->Fill(currEnergyMC, etaCurrMC);
            h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]->Fill(currEnergy, nClPerParticle);
            h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]->Fill(currEnergyMC, nClPerParticle);          
            if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
              h_clusterizer_clsspecSE_matched_E_eta[icalo][ialgo]->Fill(currEnergy, etaCurr);
              h_clusterizer_clsspecSEMC_matched_E_eta[icalo][ialgo]->Fill(currEnergyMC, etaCurrMC);     
            }
            for(int iPDG=0;iPDG<nParticlesPDG_CS;iPDG++){
              if(abs(int_CS_mcparticles_PDG[iPDG])==abs(_mcpart_PDG[mcIDCurr])){
                h_clusterizer_clsspecSE_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergy, etaCurr);
                h_clusterizer_clsspecSEMC_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergyMC, etaCurrMC);
                if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
                  h_clusterizer_clsspecSE_matched_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergy, etaCurr);
                  h_clusterizer_clsspecSEMC_matched_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergyMC, etaCurrMC);                  
                  if (trackP > 0){
                    if(h_clusterizer_clsspecSE_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]) h_clusterizer_clsspecSE_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Fill(currEnergy, currEnergy/trackP);
                    if(h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]) h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Fill(currEnergyMC, currEnergy/trackP);
                  }
                }
              }
            }
          }
        } else {
          // fill inputs for efficiency
          h_clusterizer_clsspec_E_eta[icalo][ialgo]->Fill(currEnergy, etaCurr);
          h_clusterizer_clsspecMC_E_eta[icalo][ialgo]->Fill(currEnergyMC, etaCurrMC);
          if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
            h_clusterizer_clsspec_matched_E_eta[icalo][ialgo]->Fill(currEnergy, etaCurr);
            h_clusterizer_clsspecMC_matched_E_eta[icalo][ialgo]->Fill(currEnergyMC, etaCurrMC);
          }
          if (!particleFilled && isHighest){
            h_clusterizer_clsspecSE_E_eta[icalo][ialgo]->Fill(currEnergy, etaCurr);
            h_clusterizer_clsspecSEMC_E_eta[icalo][ialgo]->Fill(currEnergyMC, etaCurrMC);          
            h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]->Fill(currEnergy, nClPerParticle);
            h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]->Fill(currEnergyMC, nClPerParticle);          
            if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
              h_clusterizer_clsspecSE_matched_E_eta[icalo][ialgo]->Fill(currEnergy, etaCurr);
              h_clusterizer_clsspecSEMC_matched_E_eta[icalo][ialgo]->Fill(currEnergyMC, etaCurrMC);          
            }
          }
          
          for(int iPDG=0;iPDG<nParticlesPDG_CS;iPDG++){
            if(abs(int_CS_mcparticles_PDG[iPDG])==abs(_mcpart_PDG[mcIDCurr])){
              h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergy, etaCurr);
              h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergyMC, etaCurrMC);
              h_clusterizer_clsspecMC_particle_eta_phi[icalo][ialgo][iPDG]->Fill(etaCurrMC,_mcpart_Phi[mcIDCurr]);
              h_clusterizer_clsspec_particle_eta_phi[icalo][ialgo][iPDG]->Fill(etaCurr, phiCurr);
              if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
                h_clusterizer_clsspec_matched_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergy, etaCurr);
                h_clusterizer_clsspecMC_matched_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergyMC, etaCurrMC);
                if (trackP > 0){
                  if(h_clusterizer_clsspec_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]) h_clusterizer_clsspec_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Fill(currEnergy, currEnergy/trackP);
                  if(h_clusterizer_clsspecMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]) h_clusterizer_clsspecMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Fill(currEnergyMC, currEnergy/trackP);
                }
              }
              if (!particleFilled && isHighest){
                h_clusterizer_clsspecSE_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergy, etaCurr);
                h_clusterizer_clsspecSEMC_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergyMC, etaCurrMC);
                if ((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched){
                  h_clusterizer_clsspecSE_matched_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergy, etaCurr);
                  h_clusterizer_clsspecSEMC_matched_particle_E_eta[icalo][ialgo][iPDG]->Fill(currEnergyMC, etaCurrMC);                  
                  if (trackP > 0){
                    if(h_clusterizer_clsspecSE_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]) h_clusterizer_clsspecSE_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Fill(currEnergy, currEnergy/trackP);
                    if(h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]) h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Fill(currEnergyMC, currEnergy/trackP);
                  }
                }
              }
            }
          }
          h_clusterizer_clsspec_PDG[icalo][ialgo]->Fill(abs(_mcpart_PDG[mcIDCurr]));
          if((_clusters_calo[ialgo][icalo].at(iclus)).cluster_isMatched) h_clusterizer_clsspec_matched_PDG[icalo][ialgo]->Fill(abs(_mcpart_PDG[mcIDCurr]));
        }
        if (!particleFilled && isHighest){
          occuranceStrct tempstructCl;
          tempstructCl.particle_ID = mcIDCurr;
          tempstructCl.highest_E = highestEnergy;
          tempstructCl.nClusters = nClPerParticle;
          uniqueIDs.push_back(tempstructCl);
//           if (tempstructCl.nClusters > 1) cout << "particle: " << tempstructCl.particle_ID << " has multiple clusters (" << tempstructCl.nClusters << ") for " << str_calorimeter[icalo].Data() << " " << str_clusterizer[ialgo].Data() << endl;
        }
      }
      uniqueIDs.clear();
      usedIDs.clear();
    }
  }
}


//*************************************************************************************************
//*************************************************************************************************
// ANCHOR save function after event loop
//*************************************************************************************************
//*************************************************************************************************
void clusterstudiesSave(){
  // gSystem->Exec("mkdir -p treeProcessing/clusterstudies");
  // define output file
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    for(int ialgo=0;ialgo<_active_algo;ialgo++){

      h_CS_clusters_NTowerMean_E[icalo][ialgo]  = new TH1D(Form("h_CS_NTowerMean_%s_%s_E",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()),"",nP, partP);
      h_clusterizer_NClMean_E[icalo][ialgo]     = new TH1D(Form("h_CS_NClMean_%s_%s_E",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()),"",nP, partP);
      h_clusterizer_NClMean_MCE[icalo][ialgo]   = new TH1D(Form("h_CS_NClMean_%s_%s_MCE",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()),"",nP, partP);

      if(h_CS_clusters_NTower_E[icalo][ialgo]){
        // make projections for JES and JER
        for (Int_t i=1; i < nP+1; i++){
          TH1D* projectionYdummy = (TH1D*)h_CS_clusters_NTower_E[icalo][ialgo]->ProjectionY(Form("projectionYdummy%d_%d_%d",icalo, ialgo,i), h_CS_clusters_NTower_E[icalo][ialgo]->GetXaxis()->FindBin(partP[i-1]+0.001), h_CS_clusters_NTower_E[icalo][ialgo]->GetXaxis()->FindBin(partP[i]+0.001),"e");
          // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
          h_CS_clusters_NTowerMean_E[icalo][ialgo]->SetBinContent(i,projectionYdummy->GetMean());
          h_CS_clusters_NTowerMean_E[icalo][ialgo]->SetBinError(i,projectionYdummy->GetStdDev());
        }
      }
      if(h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]){
        // make projections for JES and JER
        for (Int_t i=1; i < nP+1; i++){
          TH1D* projectionYdummy = (TH1D*)h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]->ProjectionY(Form("projectionYdummy%d_%d_%d",icalo, ialgo,i), h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]->GetXaxis()->FindBin(partP[i-1]+0.001), h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]->GetXaxis()->FindBin(partP[i]+0.001),"e");
          // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
          h_clusterizer_NClMean_E[icalo][ialgo]->SetBinContent(i,projectionYdummy->GetMean());
          h_clusterizer_NClMean_E[icalo][ialgo]->SetBinError(i,projectionYdummy->GetStdDev());
        }
      }
      if(h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]){
        // make projections for JES and JER
        for (Int_t i=1; i < nP+1; i++){
          TH1D* projectionYdummy = (TH1D*)h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]->ProjectionY(Form("projectionYdummy%d_%d_%d",icalo, ialgo,i), h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]->GetXaxis()->FindBin(partP[i-1]+0.001), h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]->GetXaxis()->FindBin(partP[i]+0.001),"e");
          // projectionYdummy->GetXaxis()->SetRangeUser(-0.6,0);
          h_clusterizer_NClMean_MCE[icalo][ialgo]->SetBinContent(i,projectionYdummy->GetMean());
          h_clusterizer_NClMean_MCE[icalo][ialgo]->SetBinError(i,projectionYdummy->GetStdDev());
        }
      }
    }
  }

  
  TFile* fileOutput = new TFile(Form("%s/output_CS.root",outputDir.Data()),"RECREATE");

  // write histograms
  for(int icalo=0;icalo<_active_calo;icalo++){
    if (!caloEnabled[icalo]) continue;
    Int_t caloDir           = GetCaloDirection(icalo);
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      if(h_CS_clusters_NTower_E[icalo][ialgo]) h_CS_clusters_NTower_E[icalo][ialgo]->Write();
      if(h_CS_clusters_NTowerMean_E[icalo][ialgo]) h_CS_clusters_NTowerMean_E[icalo][ialgo]->Write();
      if(h_CS_clusters_E[icalo][ialgo]) h_CS_clusters_E[icalo][ialgo]->Write();
      if(h_CS_clusters_M02_E[icalo][ialgo]) h_CS_clusters_M02_E[icalo][ialgo]->Write();
      if(h_CS_clusters_M20_E[icalo][ialgo]) h_CS_clusters_M20_E[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_matched_E_eta[icalo][ialgo]) h_clusterizer_clsspec_matched_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecMC_matched_E_eta[icalo][ialgo]) h_clusterizer_clsspecMC_matched_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_PDG[icalo][ialgo]) h_clusterizer_clsspec_PDG[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_matched_PDG[icalo][ialgo]) h_clusterizer_clsspec_matched_PDG[icalo][ialgo]->Write();
      if(h_clusterizer_1tower_E_eta[icalo][ialgo]) h_clusterizer_1tower_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_1towerMC_E_eta[icalo][ialgo]) h_clusterizer_1towerMC_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspec_E_eta[icalo][ialgo]) h_clusterizer_clsspec_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecMC_E_eta[icalo][ialgo]) h_clusterizer_clsspecMC_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecSE_E_eta[icalo][ialgo]) h_clusterizer_clsspecSE_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecSEMC_E_eta[icalo][ialgo]) h_clusterizer_clsspecSEMC_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]) h_clusterizer_clsspecSE_E_NCl[icalo][ialgo]->Write();
      if(h_clusterizer_NClMean_E[icalo][ialgo]) h_clusterizer_NClMean_E[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]) h_clusterizer_clsspecSEMC_E_NCl[icalo][ialgo]->Write();
      if(h_clusterizer_NClMean_MCE[icalo][ialgo]) h_clusterizer_NClMean_MCE[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecSE_matched_E_eta[icalo][ialgo]) h_clusterizer_clsspecSE_matched_E_eta[icalo][ialgo]->Write();
      if(h_clusterizer_clsspecSEMC_matched_E_eta[icalo][ialgo]) h_clusterizer_clsspecSEMC_matched_E_eta[icalo][ialgo]->Write();

      if(h_clusterizer_NClPerMCLabel[icalo][ialgo]) h_clusterizer_NClPerMCLabel[icalo][ialgo]->Write();
      if(h_clusterizer_ClSameMCLabel_dR[icalo][ialgo]) h_clusterizer_ClSameMCLabel_dR[icalo][ialgo]->Write();
    
      for (int iPDG = 0; iPDG < nParticlesPDG_CS; iPDG++){
        if(h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]) h_CS_clusters_M02_part_E_truth[icalo][ialgo][iPDG]->Write();
        if(h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]) h_CS_clusters_M20_part_E_truth[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][iPDG])h_clusterizer_clsspec_particle_E_eta[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][iPDG]) h_clusterizer_clsspecMC_particle_E_eta[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspecMC_particle_eta_phi[icalo][ialgo][iPDG]) h_clusterizer_clsspecMC_particle_eta_phi[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspec_particle_eta_phi[icalo][ialgo][iPDG]) h_clusterizer_clsspec_particle_eta_phi[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspecSE_particle_E_eta[icalo][ialgo][iPDG])h_clusterizer_clsspecSE_particle_E_eta[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspecSEMC_particle_E_eta[icalo][ialgo][iPDG]) h_clusterizer_clsspecSEMC_particle_E_eta[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspec_matched_particle_E_eta[icalo][ialgo][iPDG])h_clusterizer_clsspec_matched_particle_E_eta[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspecMC_matched_particle_E_eta[icalo][ialgo][iPDG]) h_clusterizer_clsspecMC_matched_particle_E_eta[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspecSE_matched_particle_E_eta[icalo][ialgo][iPDG])h_clusterizer_clsspecSE_matched_particle_E_eta[icalo][ialgo][iPDG]->Write();
        if(h_clusterizer_clsspecSEMC_matched_particle_E_eta[icalo][ialgo][iPDG]) h_clusterizer_clsspecSEMC_matched_particle_E_eta[icalo][ialgo][iPDG]->Write();
        for (Int_t et = minEtaBinCalo[caloDir]; et<maxEtaBinCalo[caloDir]; et++){
          if(h_clusterizer_clsspec_matched_particle_E_EoverP[icalo][ialgo][iPDG][et])h_clusterizer_clsspec_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Write();
          if(h_clusterizer_clsspecMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et])h_clusterizer_clsspecMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Write();
          if(h_clusterizer_clsspecSE_matched_particle_E_EoverP[icalo][ialgo][iPDG][et])h_clusterizer_clsspecSE_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Write();
          if(h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et])h_clusterizer_clsspecSEMC_matched_particle_E_EoverP[icalo][ialgo][iPDG][et]->Write();
        }
        
      }
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
