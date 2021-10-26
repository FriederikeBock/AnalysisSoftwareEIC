// ANCHOR basic struct for towers in clusterizer
typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_iL;
  int tower_trueID;
} towersStrct;

typedef struct {
  float cluster_E;
  float cluster_seed;
  float cluster_Eta;
  float cluster_Phi;
  float cluster_Z;
  float cluster_X;
  float cluster_Y;
  float cluster_M02;
  float cluster_M20;
  bool cluster_isMatched;
  vector<int> cluster_matchedTrackIDs;
  int cluster_NTowers;
  int cluster_trueID;
  int cluster_NtrueID;
} clustersStrct;

typedef struct {
  int particle_ID;
  float highest_E;
  int nClusters;
} occuranceStrct;

// enum calotype {
//   kFHCAL         = 0,
//   kFEMC         = 1,
//   kDRCALO        = 2,
//   kEEMC         = 3,
//   kCEMC         = 4,
//   kEHCAL         = 5,
//   kHCALIN       = 6,
//   kHCALOUT       = 7,
//   kLFHCAL        = 8,
//   kEEMCG         = 9,
//   kBECAL         = 10
// };
const int maxcalo = 12;
int calogeomindex[maxcalo] = {0};

TString str_calorimeter[maxcalo] = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC", "EHCAL", "HCALIN", "HCALOUT", "LFHCAL", "EEMCG", "BECAL", "FOCAL"};
int _combCalo[maxcalo]           = {kFEMC,  -1,     -1,       -1,      -1,      kEEMC,  kBECAL,   kHCALIN,    kFEMC,    -1,       -1,       -1};
int _combCalo2[maxcalo]          = {-1,     -1,     -1,       -1,      -1,      -1,     -1,       kBECAL,     -1,       -1,       -1,       -1};
int _caloTowersPhi[maxcalo]      = {0,       0,      0,        0,     100,      0,      64,       64,         0,        -1,       128,       0};
const int _active_calo = 12;

void CheckBarrelCaloVersion (){
  if (caloEnabled[kCEMC]) {
    _combCalo[kHCALIN]  = kCEMC;
    _combCalo[kHCALOUT] = kCEMC;
  }
}

enum clusterizertype {
    kMA         = 0,
    kV3         = 1,
    kV1         = 2,
    k5x5        = 3,
    kC5         = 4,
    kC3         = 5,
    k3x3        = 6,
    kDummy      = 7
};
//                      FHCAL   FEMC   DRCALO EEMC  CEMC  EHCAL HCALIN  HCALOUT LFHCAL  EEMCG BECAL
float seedE[maxcalo]  = {0.5,   0.1,   0.3,   0.1,  0.5,  0.5,  0.2,    0.5,    0.1,    0.1,  0.1   };
float aggE[maxcalo]   = {0.1,   0.005, 0.05,  0.01, 0.1,  0.1,  0.05,   0.1,    0.001,  0.01, 0.01  };

// TString str_clusterizer[7] = {"V1", "V3", "3x3", "5x5", "C3", "C5", "MA"};
TString str_clusterizer[7] = {"MA", "V3", "V1", "5x5", "C5", "C3", "3x3"};
const int maxAlgo = 7;
const int _active_algo = 1;

float _ch_DRCALO_pos_z = 1;
float _ch_FHCAL_pos_z = 1;
float _ch_FEMC_pos_z = 1;
float _ch_EHCAL_pos_z = 1;
float _ch_EEMC_pos_z = 1;
float _ch_EEMCG_pos_z = 1;
float _ch_LFHCAL_pos_z = 1;
float _ch_FOCAL_pos_z = 1;

// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }
bool acompareCl(clustersStrct lhs, clustersStrct rhs) { return lhs.cluster_E > rhs.cluster_E; }

void setINTClusterArrayToZero(int* &arrayinput){
  for(int iclus=0;iclus<_maxNclusters;iclus++){
    arrayinput[iclus]=0;
  }
}
void setFLOATClusterArrayToZero(float* &arrayinput){
  for(int iclus=0;iclus<_maxNclusters;iclus++){
    arrayinput[iclus]=0.;
  }
}
void setBOOLClusterArrayToZero(bool* &arrayinput){
  for(int iclus=0;iclus<_maxNclusters;iclus++){
    arrayinput[iclus]=false;
  }
}
// conversion functions for processing

float weightM02 = 4.5;
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, Bool_t debugOutput);

void SetGeometryIndices(){
  Long64_t nEntriesTree                 = tt_geometry->GetEntries();
  for (int ireset=0; ireset<maxcalo;ireset++) {
    calogeomindex[ireset] = -1;
  }
  for (Long64_t i=0; i<nEntriesTree;i++) {
    tt_geometry->GetEntry(i);
    calogeomindex[_calogeom_ID] = i;
  }
}

bool IsHCALCalorimeter(int caloID){
  switch (caloID){
    case kDRCALO: return true;
    case kFHCAL: return true;
    case kFEMC: return false;
    case kEHCAL: return true;
    case kEEMC: return false;
    case kHCALIN: return true;
    case kHCALOUT: return true;
    case kCEMC: return false;
    case kEEMCG: return false;
    case kLFHCAL: return true;
    case kBECAL: return false;
    case kFOCAL: return true;
    default:
      cout << "IsHCALCalorimeter: caloID " << caloID << " not defined, returning false" << endl;
      return false;
  }
  return false;
}
bool IsForwardCalorimeter(int caloID){
  switch (caloID){
    case kDRCALO: return true;
    case kFHCAL: return true;
    case kFEMC: return true;
    case kEHCAL: return true;
    case kEEMC: return true;
    case kHCALIN: return false;
    case kHCALOUT: return false;
    case kCEMC: return false;
    case kEEMCG: return true;
    case kLFHCAL: return true;
    case kBECAL: return false;
    case kFOCAL: return true;
    default:
      cout << "IsForwardCalorimeter: caloID " << caloID << " not defined, returning false" << endl;
      return false;
  }
  return false;
}

int GetCaloDirection(int caloID){
  switch (caloID){
    case kDRCALO: return 2;
    case kFHCAL: return 2;
    case kFEMC: return 2;
    case kEHCAL: return 0;
    case kEEMC: return 0;
    case kHCALIN: return 1;
    case kHCALOUT: return 1;
    case kCEMC: return 1;
    case kEEMCG: return 0;
    case kLFHCAL: return 2;
    case kBECAL: return 1;
    default:
      cout << "GetCaloDirection: caloID " << caloID << " not defined, returning -1" << endl;
      return -1;
  }
  return -1;
}


float ReturnFwdCalorimeterPosition(int caloID){
  if(IsForwardCalorimeter(caloID)){
    switch (caloID){
      case kDRCALO: return _ch_DRCALO_pos_z;
      case kFHCAL: return _ch_FHCAL_pos_z;
      case kFEMC: return _ch_FEMC_pos_z;
      case kEHCAL: return _ch_EHCAL_pos_z;
      case kEEMC: return _ch_EEMC_pos_z;
      case kEEMCG: return _ch_EEMCG_pos_z;
      case kLFHCAL: return _ch_LFHCAL_pos_z;
      case kFOCAL: return _ch_FOCAL_pos_z;
      default:
        cout << "ReturnFwdCalorimeterPosition: caloID " << caloID << " forward position not defined, returning 1" << endl;
        return 1;
    }
  }
  return 1;
}

int ReturnCaloFwdBarrelBckw(int caloID){
  // 0 for forward calos
  // 1 for barrel calos
  // 2 for backwards calos
  switch (caloID){
    case kDRCALO: return 0;
    case kFHCAL: return 0;
    case kFEMC: return 0;
    case kEHCAL: return 2;
    case kEEMC: return 2;
    case kHCALIN: return 1;
    case kHCALOUT: return 1;
    case kCEMC: return 1;
    case kEEMCG: return 2;
    case kLFHCAL: return 0;
    case kBECAL: return 1;
    case kFOCAL: return 0;
    default:
      return 1;
  }
  return 1;
}

void SetFwdCalorimeterPosition(int caloID, float zposin){
  float zpos = abs(zposin);
  switch (caloID){
    case kDRCALO: _ch_DRCALO_pos_z = zpos; break;
    case kFHCAL: _ch_FHCAL_pos_z = zpos; break;
    case kFEMC:_ch_FEMC_pos_z = zpos; break;
    case kEHCAL:_ch_EHCAL_pos_z = zpos; break;
    case kEEMC:_ch_EEMC_pos_z = zpos; break;
    case kEEMCG:_ch_EEMCG_pos_z = zpos; break;
    case kLFHCAL:_ch_LFHCAL_pos_z = zpos; break;
    case kFOCAL:_ch_FOCAL_pos_z = zpos; break;
    default: break;
  }
}

float* EtaPhiFromIndices(int ieta,int iphi,float energy = 0, int caloSelect = 0);

// ANCHOR function to return a TVector3 for the tower position based on iEta and iPhi indices
TVector3 TowerPositionVectorFromIndicesGeometry(int i_Eta,int i_Phi, int i_L, int caloSelect = 0){
  float xpos = -10000;
  float ypos = -10000;
  float zpos = -10000;

  if(calogeomindex[caloSelect]==-1) cout << "calorimeter not found in geometry!" << endl;
  tt_geometry->GetEntry(calogeomindex[caloSelect]);
  if (caloSelect == kLFHCAL){
    for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
      if(_calogeom_towers_iEta[itow]==i_Eta){
        if(_calogeom_towers_iPhi[itow]==i_Phi){
          if(_calogeom_towers_iL[itow]==i_L){
            xpos = _calogeom_towers_x[itow];
            ypos = _calogeom_towers_y[itow];
            zpos = _calogeom_towers_z[itow];
          }
        }
      }
    }
  } else {
    for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
      if(_calogeom_towers_iEta[itow]==i_Eta){
        if(_calogeom_towers_iPhi[itow]==i_Phi){
          xpos = _calogeom_towers_x[itow];
          ypos = _calogeom_towers_y[itow];
          zpos = _calogeom_towers_z[itow];
        }
      }
    }
  }
//   if(xpos==-10000 || ypos==-10000 || zpos==-10000){
//     cout << "something is terribly wrong in your geometry: x="<<xpos<<"\ty="<<ypos<<"\tz="<<zpos<<endl;
//   }
  if(IsForwardCalorimeter(caloSelect) && zpos!=-10000){
    if(ReturnFwdCalorimeterPosition(caloSelect)==1){
      SetFwdCalorimeterPosition(caloSelect, zpos);
    }
  }
  TVector3 twrPositionVec(xpos,ypos,zpos);
  return twrPositionVec;
}

int ReturnMaxTowerCalo(int caloID){
  switch (caloID){
    case kDRCALO: return _nTowers_DRCALO;
    case kFHCAL: return _nTowers_FHCAL;
    case kFEMC: return _nTowers_FEMC;
    case kEHCAL: return _nTowers_EHCAL;
    case kEEMC: return _nTowers_EEMC;
    case kHCALIN: return _nTowers_HCALIN;
    case kHCALOUT: return _nTowers_HCALOUT;
    case kCEMC: return _nTowers_CEMC;
    case kEEMCG: return _nTowers_EEMCG;
    case kLFHCAL: return _nTowers_LFHCAL;
    case kBECAL: return _nTowers_BECAL;
    case kFOCAL: return _nTowers_FOCAL;
    default:
      cout << "ReturnMaxTowerCalo: caloID " << caloID << " not defined, returning -1" << endl;
      return -1;
  }
  return -1;
}



int ReturnProjectionIndexForCalorimeter(int caloID,  bool alternate){
  if (!alternate){
    switch (caloID){
      case kDRCALO: return 1;
      case kFHCAL: return 5;
      case kFEMC: return 6;
      case kEHCAL: return 60;
      case kEEMC: return 61;
      case kHCALIN: return 62;
      case kHCALOUT: return 63;
      case kCEMC: return 64;
      case kEEMCG: return 65;
      case kLFHCAL: return 67;
      case kBECAL: return 66;
      case kFOCAL: return 85;
      default:
        cout << "ReturnProjectionIndexForCalorimeter: caloID " << caloID << " not defined, returning -1" << endl;
        return -1;
    }
  } else {
    switch (caloID){
      case kDRCALO: return 1;
      case kFHCAL: return 1;
      case kFEMC: return 1;
      case kEHCAL: return 4;
      case kEEMC: return 4;
      case kHCALIN: return 8;
      case kHCALOUT: return 8;
      case kCEMC: return 8;
      case kEEMCG: return 4;
      case kLFHCAL: return 1;
      case kBECAL: return 8;
      default:
        cout << "ReturnProjectionIndexForCalorimeter: caloID " << caloID << " not defined, returning -1" << endl;
        return -1;
    }
    
  }
  return -1;
}

int ReturnCalorimeterFromProjectionIndex(int projID){
  switch (projID){
    case 1: return kDRCALO;
    case 5: return kFHCAL;
    case 6: return kFEMC;
    case 60: return kEHCAL;
    case 61: return kEEMC;
    case 62: return kHCALIN;
    case 63: return kHCALOUT;
    case 64: return kCEMC;
    case 65: return kEEMCG;
    case 67: return kLFHCAL;
    case 66: return kBECAL;
    case 85: return kFOCAL;
    default:
      // cout << "ReturnCalorimeterFromProjectionIndex: projID " << projID << " not defined, returning -1" << endl;
      return -1;
  }
  return -1;
}

float ReturnTrackMatchingWindowForCalo(int caloID){
  switch (caloID){
    case kDRCALO: return 5;
    case kFHCAL: return 10;
    case kFEMC: return 5.5;
    case kEHCAL: return 10;
    case kEEMC: return 4.0;
    case kHCALIN: return 0.1;
    case kHCALOUT: return 0.1;
    case kCEMC: return 0.05;
    case kEEMCG: return 7;
    case kLFHCAL: return 10;
    case kBECAL: return 0.2;
    case kFOCAL: return 3;
    default:
      cout << "ReturnTrackMatchingWindowForCalo: caloID " << caloID << " not defined, returning -1" << endl;
      return -1;
  }
  return -1;
}


// ANCHOR function to determine shower shape
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, bool debugOutput){
    static float returnVariables[7]; //0:M02, 1:M20, 2:eta, 3: phi
    float w_tot = 0;
    std::vector<float> w_i;
    TVector3 vecTwr;
    TVector3 vecTwrTmp;
    float zHC = 1;

    vecTwr = {0.,0.,0.};
    //calculation of weights and weighted position vector
    for(int cellI=0; cellI<(int)cluster_towers.size(); cellI++){
        w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
        w_tot += w_i.at(cellI);
        vecTwrTmp = TowerPositionVectorFromIndicesGeometry(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, cluster_towers.at(cellI).tower_iL, caloSelect);
//         cout << caloSelect << "\t" << vecTwrTmp.X() << "\t" << vecTwrTmp.Y() << "\t" << vecTwrTmp.Z() << endl;
        vecTwr += w_i.at(cellI)*vecTwrTmp;
        if(cellI==0 && ReturnFwdCalorimeterPosition(caloSelect))zHC=vecTwrTmp.Z();
    }
    returnVariables[2]=vecTwr.Eta();
    returnVariables[3]=vecTwr.Phi(); //(vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());
//     cout << "X: "<< vecTwr.X() << "\t" << " Y: "<< vecTwr.Y() << "\t" << " Z: "<< vecTwr.Z() << "\t zHC: " <<  zHC << endl;
    if(IsForwardCalorimeter(caloSelect)){
      vecTwr*=abs(zHC/vecTwr.Z()); // wtot/Ntowers?
    }
    returnVariables[4]=vecTwr.X();
    returnVariables[5]=vecTwr.Y();
    returnVariables[6]=vecTwr.Z();

    //calculation of M02
    float delta_phi_phi[4] = {0};
    float delta_eta_eta[4] = {0};
    float delta_eta_phi[4] = {0};
    int supModLeadingCell = -1;
    for(int cellI=0; cellI<(int)cluster_towers.size(); cellI++){
        int iphi=cluster_towers.at(cellI).tower_iPhi;
        int ieta=cluster_towers.at(cellI).tower_iEta;
        delta_phi_phi[1] += (w_i.at(cellI)*iphi*iphi)/w_tot;
        delta_phi_phi[2] += (w_i.at(cellI)*iphi)/w_tot;
        delta_phi_phi[3] += (w_i.at(cellI)*iphi)/w_tot;

        delta_eta_eta[1] += (w_i.at(cellI)*ieta*ieta)/w_tot;
        delta_eta_eta[2] += (w_i.at(cellI)*ieta)/w_tot;
        delta_eta_eta[3] += (w_i.at(cellI)*ieta)/w_tot;

        delta_eta_phi[1] += (w_i.at(cellI)*ieta*iphi)/w_tot;
        delta_eta_phi[2] += (w_i.at(cellI)*iphi)/w_tot;
        delta_eta_phi[3] += (w_i.at(cellI)*ieta)/w_tot;
    }
    delta_phi_phi[0] = delta_phi_phi[1] - (delta_phi_phi[2] * delta_phi_phi[3]);
    delta_eta_eta[0] = delta_eta_eta[1] - (delta_eta_eta[2] * delta_eta_eta[3]);
    delta_eta_phi[0] = delta_eta_phi[1] - (delta_eta_phi[2] * delta_eta_phi[3]);

    float calcM02 = 0.5 * ( delta_phi_phi[0] + delta_eta_eta[0] ) + TMath::Sqrt( 0.25 * TMath::Power( ( delta_phi_phi[0] - delta_eta_eta[0] ), 2 ) + TMath::Power( delta_eta_phi[0], 2 ) );
    float calcM20 = 0.5 * ( delta_phi_phi[0] + delta_eta_eta[0] ) - TMath::Sqrt( 0.25 * TMath::Power( ( delta_phi_phi[0] - delta_eta_eta[0] ), 2 ) + TMath::Power( delta_eta_phi[0], 2 ) );
    if(debugOutput) cout << "M02_calc: " << calcM02 << "\t\t = 0.5 * ( " << delta_phi_phi[0] <<" + "<<delta_eta_eta[0]<<" ) + TMath::Sqrt( 0.25 * TMath::Power( ( "<<delta_phi_phi[0]<<" - "<<delta_eta_eta[0]<<" ), 2 ) + TMath::Power( "<<delta_eta_phi[0]<<", 2 ) ) "<< endl;
    returnVariables[0]=calcM02;
    returnVariables[1]=calcM20;
    return returnVariables;

}


// 3x3 global cluster variables
std::vector<clustersStrct> _clusters_calo[maxAlgo][maxcalo];

bool loadClusterizerInput(
  int clusterizerEnum,
  int caloEnum
){
  // cout << clusterizerEnum << "\t" << caloEnum << endl;
  if (clusterizerEnum != kV1 ){
    if (!_clusters_calo[clusterizerEnum][caloEnum].empty() && caloEnabled[caloEnum])
      return true;
    else 
      return false;
  } else {
    
    if(caloEnum==kFHCAL){      
      if (!_clusters_calo[clusterizerEnum][caloEnum].empty() && caloEnabled[caloEnum]){
        return true;
      } else {
        if (_nclusters_FHCAL > 0 ){
          for (Int_t iclus = 0; iclus < (Int_t)_nclusters_FHCAL; iclus++){
            clustersStrct tempstructC;
            tempstructC.cluster_seed    = -1;
            tempstructC.cluster_E       =_clusters_FHCAL_E[iclus];
            tempstructC.cluster_NTowers =_clusters_FHCAL_NTower[iclus];
            tempstructC.cluster_M02     = 0;
            tempstructC.cluster_M20     = 0;
            tempstructC.cluster_Eta     = _clusters_FHCAL_Eta[iclus];
            tempstructC.cluster_Phi     = _clusters_FHCAL_Phi[iclus];
            tempstructC.cluster_X       = 0;
            tempstructC.cluster_Y       = 0;
            tempstructC.cluster_Z       = 0;
            tempstructC.cluster_isMatched = false;
            tempstructC.cluster_trueID  = _clusters_FHCAL_trueID[iclus];
            tempstructC.cluster_NtrueID = 0;
            _clusters_calo[clusterizerEnum][caloEnum].push_back(tempstructC);
          } 
          return true;
        } else {
          return false;
        }
      }
    } else if(caloEnum==kFEMC){
      if (!_clusters_calo[clusterizerEnum][caloEnum].empty() && caloEnabled[caloEnum]){
        return true;
      } else {
        if (_nclusters_FEMC > 0 ){
          for (Int_t iclus = 0; iclus < (Int_t)_nclusters_FEMC; iclus++){
            clustersStrct tempstructC;
            tempstructC.cluster_seed    = -1;
            tempstructC.cluster_E       =_clusters_FEMC_E[iclus];
            tempstructC.cluster_NTowers =_clusters_FEMC_NTower[iclus];
            tempstructC.cluster_M02     = 0;
            tempstructC.cluster_M20     = 0;
            tempstructC.cluster_Eta     = _clusters_FEMC_Eta[iclus];
            tempstructC.cluster_Phi     = _clusters_FEMC_Phi[iclus];
            tempstructC.cluster_X       = 0;
            tempstructC.cluster_Y       = 0;
            tempstructC.cluster_Z       = 0;
            tempstructC.cluster_isMatched = false;
            tempstructC.cluster_trueID  = _clusters_FEMC_trueID[iclus];
            tempstructC.cluster_NtrueID = 0;
            _clusters_calo[clusterizerEnum][caloEnum].push_back(tempstructC);
          } 
          return true;
        } else {
          return false;
        }
      } 
    } else {
      if (!_clusters_calo[clusterizerEnum][caloEnum].empty() && caloEnabled[caloEnum])
        return true;
      else 
        return false;
    }
  }

}

float getCalibrationValue(float clusterE, int caloEnum, int algoEnum){
  if(caloEnum == kFHCAL){
    float addCorr = 1.03;
    if(algoEnum==kV1){
      float paramscalib[5]= { 2.29669, 0.48742, 2.09749, -703.543, 469.853 }; //noshift
      return (addCorr * ( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kV3){
      float paramscalib[5]= { 4.54775, 1.59832, 6.90497, -10110.7, 11986.3 }; //noshift
      return (addCorr * ( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==k3x3){
      float paramscalib[5]= { 0.403972, 0.324973, 2.77029, -170466, -578827 }; //noshift
      return (addCorr * ( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==k5x5){
      float paramscalib[5]= { 31.7632, 5.39668, 78.4091, -83027.6, 556743 }; //noshift
      return (addCorr * ( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kC3){
      float paramscalib[5]= { -3.77545, 4.13685, 2.15795, -656.349, 301.55 }; //noshift
      return (addCorr * ( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kC5){
      float paramscalib[5]= { 97.7482, 17.6935, 262.447, -191344, 1.51638e+06 }; //noshift
      return (addCorr * ( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kMA){
      float paramscalib[5]= { 3.02343, 0.408232, 2.49633, -2426.69, 3259.87 }; //noshift
      return (addCorr * ( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else {
      return 1.0;
    }
  } else if(caloEnum == kFEMC){
    if(algoEnum==kV1){
      float paramscalib[5]= {3.61545, 0.0317621, 1.80493, -1159.32, 1869.31}; //noshift
      return (( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kV3){
      float paramscalib[5]= {3.48283, 0.0507623, 1.94421, -1062.14, 1689.39}; //noshift
      return (( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==k3x3){
      float paramscalib[5]= {3.66257, 0.0291583, 2.09469, -1799.62, 2795.61}; //noshift
      return (( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==k5x5){
      float paramscalib[5]= {3.64167, 0.0288047, 2.06243, -2506.92, 3843.72}; //noshift
      return (( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kC3){
      float paramscalib[5]= {3.85225, 0.0153552, 2.22351, -2221.65, 3445.8}; //noshift
      return (( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kC5){
      float paramscalib[5]= {3.78109, 0.0301875, 2.18153, -2171.59, 3377.11}; //noshift
      return (( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else if(algoEnum==kMA){
      float paramscalib[5]= {3.43486, 0.120574, 1.81206, -613.209, 913.285}; //noshift
      return (( paramscalib[0] + paramscalib[1] * TMath::Log(clusterE) ) / ( 1 + ( paramscalib[2] * TMath::Exp( ( clusterE - paramscalib[3] ) / paramscalib[4] ) ) ));
    } else {
      return 1.0;
    }
  } else {
    return 1.0;
  }
  return 1.0;
}

float getEnergySmearing( int caloEnum, int algoEnum){
  _fRandom.SetSeed(0);
  if(caloEnum==kFHCAL){
    if(algoEnum==kMA){
      return _fRandom.Gaus(1,0.0985);
    } else if(algoEnum==kV1){
      return _fRandom.Gaus(1,0.10);
    } else if(algoEnum==kV3){
      return _fRandom.Gaus(1,0.04);
    } else {
      return 1.0;
    }
  } else {
    return 1.0;
  }
  return 1.0;
}

void fillHCalClustersIntoJetFindingInputs( int caloEnum, int clusterizerEnum,  
  std::vector<float> & jetf_hcal_E, std::vector<float> & jetf_hcal_px, std::vector<float> & jetf_hcal_py, std::vector<float> & jetf_hcal_pz,
  std::vector<float> & jetf_calo_E, std::vector<float> & jetf_calo_px, std::vector<float> & jetf_calo_py, std::vector<float> & jetf_calo_pz,
  std::vector<float> & jetf_all_E, std::vector<float> & jetf_all_px, std::vector<float> & jetf_all_py, std::vector<float> & jetf_all_pz
)
{
  for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[clusterizerEnum][caloEnum].size(); iclus++){
      if(!(_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_isMatched){
          double pt = (_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E / cosh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);
          double px = pt * cos((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double py = pt * sin((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double pz = pt * sinh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);
          jetf_hcal_px.push_back(px);
          jetf_hcal_py.push_back(py);
          jetf_hcal_pz.push_back(pz);
          jetf_hcal_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_calo_px.push_back(px);
          jetf_calo_py.push_back(py);
          jetf_calo_pz.push_back(pz);
          jetf_calo_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_all_px.push_back(px);
          jetf_all_py.push_back(py);
          jetf_all_pz.push_back(pz);
          jetf_all_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
      }
  }
}

void fillECalClustersIntoJetFindingInputs(
  int caloEnum, int clusterizerEnum,  
  std::vector<float> & jetf_emcal_E, std::vector<float> & jetf_emcal_px, std::vector<float> & jetf_emcal_py, std::vector<float> & jetf_emcal_pz,
  std::vector<float> & jetf_calo_E, std::vector<float> & jetf_calo_px, std::vector<float> & jetf_calo_py, std::vector<float> & jetf_calo_pz,
  std::vector<float> & jetf_all_E, std::vector<float> & jetf_all_px, std::vector<float> & jetf_all_py, std::vector<float> & jetf_all_pz,
  std::vector<float> & jetf_full_E, std::vector<float> & jetf_full_px, std::vector<float> & jetf_full_py, std::vector<float> & jetf_full_pz
)
{
    for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[clusterizerEnum][caloEnum].size(); iclus++){
      if(!(_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_isMatched){
          double pt = (_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E / cosh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);
          double px = pt * cos((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double py = pt * sin((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Phi);
          double pz = pt * sinh((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_Eta);
          jetf_emcal_px.push_back(px);
          jetf_emcal_py.push_back(py);
          jetf_emcal_pz.push_back(pz);
          jetf_emcal_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_calo_px.push_back(px);
          jetf_calo_py.push_back(py);
          jetf_calo_pz.push_back(pz);
          jetf_calo_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_all_px.push_back(px);
          jetf_all_py.push_back(py);
          jetf_all_pz.push_back(pz);
          jetf_all_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
          jetf_full_px.push_back(px);
          jetf_full_py.push_back(py);
          jetf_full_pz.push_back(pz);
          jetf_full_E.push_back((_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_E);
      }
  }
}
