#ifndef CALOHEADER_H
#define CALOHEADER_H

// ANCHOR basic struct for towers in clusterizer
struct towersStrct{
  towersStrct(): tower_E(0), tower_iEta(-1), tower_iPhi(-1), tower_iL(-1), tower_trueID(-10000) {}
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_iL;
  int tower_trueID;
} ;

struct matchingStrct {
  matchingStrct(): id(-1), dPhi(-1000), dEta(-1000), dX(-1000), dY(-1000) {}
  int id;
  float dPhi;
  float dEta;
  float dX;
  float dY;
};

struct clustersStrct{
  clustersStrct(): cluster_E(0.), cluster_seed(0.), cluster_Eta(-10.), cluster_Phi(-10.), cluster_X(0.) , cluster_Y(0.), cluster_Z(0.), cluster_M02(0.), cluster_M20(0.), cluster_isMatched(false), cluster_isMatchedECal(false), cluster_isMatchedHCal(false), cluster_NTowers(0), cluster_trueID(-10000), cluster_NtrueID(0) {}
  float cluster_E;
  float cluster_seed;
  float cluster_Eta;
  float cluster_Phi;
  float cluster_X;
  float cluster_Y;
  float cluster_Z;
  float cluster_M02;
  float cluster_M20;
  bool cluster_isMatched;
  bool cluster_isMatchedECal;
  bool cluster_isMatchedHCal;
  int cluster_NTowers;
  int cluster_trueID;
  int cluster_NtrueID;
  std::vector<matchingStrct> cluster_matchedTracks;
  std::vector<matchingCalStrct> cluster_matchedECals;
  std::vector<matchingCalStrct> cluster_matchedHCals;
  std::vector<towersStrct> cluster_towers;
} ;

struct occuranceStrct{
  occuranceStrct(): particle_ID(-10000), highest_E(0.), nClusters(0) {}
  int particle_ID;
  float highest_E;
  int nClusters;
} ;

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
int _combCalo[maxcalo]           = {kFEMC,  -1,     -1,       -1,      -1,      kEEMC,  kBECAL,   kBECAL,     kFEMC,    -1,       -1,       -1};
int _combCalo2[maxcalo]          = {-1,     -1,     -1,       -1,      -1,      -1,     -1,       -1,         -1,       -1,       -1,       -1};
// int _combCalo2[maxcalo]          = {-1,     -1,     -1,       -1,      -1,      -1,     -1,       kHCALIN,    -1,       -1,       -1,       -1};
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
// TString str_clusterizer[7] = {"V3", "V3", "V1", "5x5", "C5", "C3", "3x3"};
const int maxAlgo = 7;
const int _active_algo = 1;

enum calibBECAL{
    kSciGlas      = 0,
    kPbGlasG4     = 1,
    kPbGlasTF1    = 2
};

calibBECAL settingCalibBECAL = kSciGlas;

float _ch_DRCALO_pos_z = 1;
float _ch_FHCAL_pos_z = 1;
float _ch_FEMC_pos_z = 1;
float _ch_EHCAL_pos_z = 1;
float _ch_EEMC_pos_z = 1;
float _ch_EEMCG_pos_z = 1;
float _ch_LFHCAL_pos_z = 1;
float _ch_FOCAL_pos_z = 1;

TGraphErrors* graphsM02Cuts[maxcalo]  = {NULL};

// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }
bool acompareTrueID(towersStrct lhs, towersStrct rhs) { return lhs.tower_trueID > rhs.tower_trueID; }
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
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float cluster_E_calc, int caloSelect, Bool_t debugOutput);


const int _maxNtowers1D = 200;
const int _maxNtowersL = 10;
TVector3 caloPositionArray[maxcalo][_maxNtowers1D/*iEta*/][_maxNtowers1D/*iPhi*/][_maxNtowersL/*iL*/];
void SetGeometryIndices(){
  Long64_t nEntriesTree                 = tt_geometry->GetEntries();
  for (int ireset=0; ireset<maxcalo;ireset++) {
    calogeomindex[ireset] = -1;
  }
  for (Long64_t i=0; i<nEntriesTree;i++) {
    tt_geometry->GetEntry(i);
    calogeomindex[_calogeom_ID] = i;

    for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
      caloPositionArray[_calogeom_ID][_calogeom_towers_iEta[itow]][_calogeom_towers_iPhi[itow]][_calogeom_towers_iL[itow]] = TVector3(_calogeom_towers_x[itow],_calogeom_towers_y[itow],_calogeom_towers_z[itow]);
    }
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
      std::cout << "IsHCALCalorimeter: caloID " << caloID << " not defined, returning false" << std::endl;
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
      std::cout << "IsForwardCalorimeter: caloID " << caloID << " not defined, returning false" << std::endl;
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
      std::cout << "GetCaloDirection: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}

float GetWeightCaloShowerShape(int caloID){
  // increasing the weight makes the distributions wider
  switch (caloID){
    case kDRCALO: return 4.5;
    case kFHCAL: return 4.5;
    case kFEMC: return 3.5;
    case kEHCAL: return 4.5;
    case kEEMC: return 4.5;
    case kHCALIN: return 4.5;
    case kHCALOUT: return 4.5;
    case kCEMC: return 4.5;
    case kEEMCG: return 4.5;
    case kLFHCAL: return 4.5;
    case kBECAL: return 4.;
    default:
      std::cout << "GetCaloDirection: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}

void LoadM02CutValueGraphs(TString nameFile){
  TFile* fileTemp = new TFile(nameFile.Data(), "READ");
  for (Int_t icalo = 0; icalo < maxcalo; icalo++){
    graphsM02Cuts[icalo] = (TGraphErrors*)fileTemp->Get(Form("%s/M02CutVal90_%s", str_calorimeter[icalo].Data(), str_calorimeter[icalo].Data() ));
    if (graphsM02Cuts[icalo]){
      std::cout << "found M02 cut graph for "<< str_calorimeter[icalo].Data()<< endl;
      graphsM02Cuts[icalo]->Print();
    }
  }
}

float GetM02CutValueForEMC(int caloID, float clusterE){
  switch (caloID){
    case kFEMC: 
      if (graphsM02Cuts[caloID]){
        if (clusterE > graphsM02Cuts[caloID]->GetX()[graphsM02Cuts[caloID]->GetN()-1])
          return graphsM02Cuts[caloID]->GetY()[graphsM02Cuts[caloID]->GetN()-1];
        else if (clusterE < graphsM02Cuts[caloID]->GetX()[0])
          return graphsM02Cuts[caloID]->GetY()[0];
        else 
          return graphsM02Cuts[caloID]->Eval(clusterE);
      } else {
        return 1;
      }
      break;
    case kEEMC: 
      if (graphsM02Cuts[caloID]){
        if (clusterE > graphsM02Cuts[caloID]->GetX()[graphsM02Cuts[caloID]->GetN()-1])
          return graphsM02Cuts[caloID]->GetY()[graphsM02Cuts[caloID]->GetN()-1];
        else if (clusterE < graphsM02Cuts[caloID]->GetX()[0])
          return graphsM02Cuts[caloID]->GetY()[0];
        else 
          return graphsM02Cuts[caloID]->Eval(clusterE);
      } else {
        return 1;
      }
      break;
    case kEEMCG: return 1;
    case kBECAL: 
      if (graphsM02Cuts[caloID]){
        if (clusterE > graphsM02Cuts[caloID]->GetX()[graphsM02Cuts[caloID]->GetN()-1])
          return graphsM02Cuts[caloID]->GetY()[graphsM02Cuts[caloID]->GetN()-1];
        else if (clusterE < graphsM02Cuts[caloID]->GetX()[0])
          return graphsM02Cuts[caloID]->GetY()[0];
        else 
          return graphsM02Cuts[caloID]->Eval(clusterE);
      } else {
        return 1;
      }
      break;
    default:
      std::cout << "GetM02CutValueForEMC: caloID " << caloID << " not defined, returning 1000" << std::endl;
      return 1000;
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
        std::cout << "ReturnFwdCalorimeterPosition: caloID " << caloID << " forward position not defined, returning 1" << std::endl;
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

float CorrectEtaPositionCalo (int caloID, float etarec){
  switch(caloID){
    case kDRCALO: 
    case kFHCAL: 
    case kFEMC:
    case kEEMC:
    case kEEMCG:
    case kFOCAL:  
    case kHCALIN: 
    case kHCALOUT:
    case kBECAL: 
    case kCEMC: 
      return etarec;
    case kEHCAL:
      return etarec+0.03;
    case kLFHCAL:
      return etarec-0.04;
    default:
      return etarec;
  }
}


float* EtaPhiFromIndices(int ieta,int iphi,float energy = 0, int caloSelect = 0);


// ANCHOR function to return a TVector3 for the tower position based on iEta and iPhi indices
TVector3 TowerPositionVectorFromIndicesGeometry(int i_Eta,int i_Phi, int i_L, int caloSelect = 0){
  TVector3 twrPositionVec = caloPositionArray[caloSelect][i_Eta][i_Phi][caloSelect == kLFHCAL ? i_L : -1];
  // cout << caloSelect << "\t" << i_Eta<< "\t" << i_Phi<< "\t" << i_L << endl;
  // cout << "\tvec: " << twrPositionVec.x() << "\t" << twrPositionVec.y()<< "\t" << twrPositionVec.z() << endl;
  if(IsForwardCalorimeter(caloSelect) && twrPositionVec.z()!=-10000){
    if(ReturnFwdCalorimeterPosition(caloSelect)==1){
      SetFwdCalorimeterPosition(caloSelect, twrPositionVec.z());
    }
  }
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
      std::cout << "ReturnMaxTowerCalo: caloID " << caloID << " not defined, returning -1" << std::endl;
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
      case kLFHCAL: return 5;
      case kBECAL: return 66;
      case kFOCAL: return 85;
      default:
        std::cout << "ReturnProjectionIndexForCalorimeter: caloID " << caloID << " not defined, returning -1" << std::endl;
        return -1;
    }
  } else {
    switch (caloID){
      case kDRCALO: return 1;
      case kFHCAL: return 1;
      case kFEMC: return 1;
      case kEHCAL: return 4;
      case kEEMC: return 4;
      case kHCALIN: return 7;
      case kHCALOUT: return 7;
      case kCEMC: return 7;
      case kEEMCG: return 4;
      case kLFHCAL: return 1;
      case kBECAL: return 7;
      default:
        std::cout << "ReturnProjectionIndexForCalorimeter: caloID " << caloID << " not defined, returning -1" << std::endl;
        return -1;
    }
    
  }
  return -1;
}

// ***********************************************************************************
// matching windows for tracks to calo's going outward dX/dy for forward, dPhi/dEta for barrel
// ***********************************************************************************
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
    case kBECAL: return 0.05;
    case kFOCAL: return 3;
    default:
      std::cout << "ReturnTrackMatchingWindowForCalo: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}

// ***********************************************************************************
// matching windows for Calo's going inward in dPhi
// ***********************************************************************************
float ReturnPhiCaloMatching(int caloID, int caloID2){
  switch (caloID){
    case kDRCALO: return 0.1;
    case kFHCAL: return 0.1;
    case kFEMC: return 0.1;
    case kEHCAL: return 0.1;
    case kEEMC: return 0.1;
    case kHCALIN: return 0.1;
    case kHCALOUT: 
      switch (caloID2){
        case kBECAL:
          return 0.08;
        case kHCALIN:
          return 0.02;
        default:
          return 0.1;
      }
    case kCEMC: return 0.1;
    case kEEMCG: return 0.1;
    case kLFHCAL: return 0.07;
    case kBECAL: return 0.05;
    case kFOCAL: return 0.1;
    default:
      std::cout << "ReturnPhiCaloMatching: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}

// ***********************************************************************************
// matching windows for Calo's going inward in dEta
// ***********************************************************************************
float ReturnEtaCaloMatching(int caloID, int caloID2){
  switch (caloID){
    case kDRCALO: return 0.1;
    case kFHCAL: return 0.1;
    case kFEMC: return 0.1;
    case kEHCAL: return 0.1;
    case kEEMC: return 0.1;
    case kHCALIN: return 0.1;
    case kHCALOUT: 
      switch (caloID2){
        case kBECAL:
          return 0.05;
        case kHCALIN:
          return 0.02;
        default:
          return 0.1;
      }
    case kCEMC: return 0.1;
    case kEEMCG: return 0.1;
    case kLFHCAL: return 0.1;
    case kBECAL: return 0.05;
    case kFOCAL: return 0.1;
    default:
      std::cout << "ReturnEtaCaloMatching: caloID " << caloID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}


// ANCHOR function to determine shower shape
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float cluster_E_calc, int caloSelect, bool debugOutput){
    static float returnVariables[8]; //0:M02, 1:M20, 2:eta, 3: phi
    float w_tot = 0;
    std::vector<float> w_i;
    TVector3 vecTwr;
    TVector3 vecTwrTmp;
    float zHC     = 1;
    float w_0     = GetWeightCaloShowerShape(caloSelect);
    
    vecTwr = {0.,0.,0.};
    //calculation of weights and weighted position vector
    for(int cellI=0; cellI<(int)cluster_towers.size(); cellI++){
        w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
        w_tot += w_i.at(cellI);
        vecTwrTmp = TowerPositionVectorFromIndicesGeometry(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, cluster_towers.at(cellI).tower_iL, caloSelect);
//         std::cout << caloSelect << "\t" << vecTwrTmp.X() << "\t" << vecTwrTmp.Y() << "\t" << vecTwrTmp.Z() << std::endl;
        vecTwr += w_i.at(cellI)*vecTwrTmp;
        if(cellI==0 && ReturnFwdCalorimeterPosition(caloSelect))zHC=vecTwrTmp.Z();
    }
    // correct Eta position for average shift in calo 
    returnVariables[2]= CorrectEtaPositionCalo(caloSelect, vecTwr.Eta());
    returnVariables[3]= vecTwr.Phi(); //(vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());
//     std::cout << "X: "<< vecTwr.X() << "\t" << " Y: "<< vecTwr.Y() << "\t" << " Z: "<< vecTwr.Z() << "\t zHC: " <<  zHC << std::endl;
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
    float dispersion = 0;
    int supModLeadingCell = -1;
    for(int cellI=0; cellI<(int)cluster_towers.size(); cellI++){
//       vecTwrTmp = TowerPositionVectorFromIndicesGeometry(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, cluster_towers.at(cellI).tower_iL, caloSelect);
//       int iphi = vecTwrTmp.Phi();
//       int ieta = vecTwrTmp.Eta();      
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

      vecTwrTmp = TowerPositionVectorFromIndicesGeometry(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, cluster_towers.at(cellI).tower_iL, caloSelect);
      if(IsForwardCalorimeter(caloSelect)){
        // scale cluster position to z-plane
        vecTwr*=abs(vecTwrTmp.Z()/vecTwr.Z());
      } else {
        // scale cluster position to same radial position as towers
        vecTwr*=abs(sqrt(pow(vecTwrTmp.X(),2)+pow(vecTwrTmp.Y(),2))/sqrt(pow(vecTwr.X(),2)+pow(vecTwr.Y(),2)));
      }
      float dx2 = pow(vecTwrTmp.X()-vecTwr.X(),2);
      float dy2 = pow(vecTwrTmp.Y()-vecTwr.Y(),2);
      float dz2 = pow(vecTwrTmp.Z()-vecTwr.Z(),2);
      dispersion+= (w_i.at(cellI)*(dx2+dy2+dz2))/w_tot;
    }
    returnVariables[7]=dispersion;
    delta_phi_phi[0] = delta_phi_phi[1] - (delta_phi_phi[2] * delta_phi_phi[3]);
    delta_eta_eta[0] = delta_eta_eta[1] - (delta_eta_eta[2] * delta_eta_eta[3]);
    delta_eta_phi[0] = delta_eta_phi[1] - (delta_eta_phi[2] * delta_eta_phi[3]);

    float calcM02 = 0.5 * ( delta_phi_phi[0] + delta_eta_eta[0] ) + TMath::Sqrt( 0.25 * TMath::Power( ( delta_phi_phi[0] - delta_eta_eta[0] ), 2 ) + TMath::Power( delta_eta_phi[0], 2 ) );
    float calcM20 = 0.5 * ( delta_phi_phi[0] + delta_eta_eta[0] ) - TMath::Sqrt( 0.25 * TMath::Power( ( delta_phi_phi[0] - delta_eta_eta[0] ), 2 ) + TMath::Power( delta_eta_phi[0], 2 ) );
    if(debugOutput) std::cout << "M02_calc: " << calcM02 << "\t\t = 0.5 * ( " << delta_phi_phi[0] <<" + "<<delta_eta_eta[0]<<" ) + TMath::Sqrt( 0.25 * TMath::Power( ( "<<delta_phi_phi[0]<<" - "<<delta_eta_eta[0]<<" ), 2 ) + TMath::Power( "<<delta_eta_phi[0]<<", 2 ) ) "<< std::endl;
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
  // std::cout << clusterizerEnum << "\t" << caloEnum << std::std::endl;
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

float getCalibrationValue(float clusterE, int caloEnum, int algoEnum, int MCtrueID){
  if(caloEnum == kFHCAL){
      return 1.0;
  } else if(caloEnum == kEEMC){
    if(algoEnum==kMA){
//       return ( 2.21528e+00 + 7.34167e-02 * TMath::Log(clusterE) ) / ( 1 + ( 3.41200e-01 * TMath::Exp( ( clusterE + 4.66905e+02 ) / 3.10970e+02 ) ) ); //2mm Carbon
      double addCalibFromPi0 = 1;
      return ( 2.21528e+00 + 7.34167e-02 * TMath::Log(clusterE) ) / ( 1 + ( 3.41200e-01 * TMath::Exp( ( clusterE + 4.66905e+02 ) / 3.10970e+02 ) ) )*1.04*addCalibFromPi0;  // 0.5mm Carbon
    } else {
      return 0.9;
    }
    return 1;
  } else if(caloEnum == kBECAL){
    if(algoEnum==kMA){
      if ( settingCalibBECAL == kSciGlas){
        float corrFac1  = ( 5.91169e+06 + 1.54628e+05 * TMath::Log(clusterE) ) / ( 1 + ( 2.14932e+06 * TMath::Exp( ( clusterE + 1.22185e+02 ) / 1.09576e+02 ) ) );
        return corrFac1;        
      } else if ( settingCalibBECAL == kPbGlasG4){
        float params[5] = {1.1, 0.0161174, 0.219832, 3.29418, 27.857 };
        float corrFac1  = ( params[0] + params[1] * TMath::Log(clusterE) ) / ( 1 + ( params[2] * TMath::Exp( ( clusterE  + params[3] ) / params[4] ) ) );
        return corrFac1*1.045;
      } else if ( settingCalibBECAL == kPbGlasTF1){
        float params[5] = {1.57990e+07, 5.17276e+05, 4.88518e+06, 1.39930e+02, 1.08064e+02};
        float corrFac1  = ( params[0] + params[1] * TMath::Log(clusterE) ) / ( 1 + ( params[2] * TMath::Exp( ( clusterE  + params[3] ) / params[4] ) ) );
        return corrFac1;
      }
    } else {
      return 1.;
    }
    return 1;
  } else if(caloEnum == kFEMC){
      double addCalibFromPi0 = 0.12/0.135;
      return ( 1.85096e+00  + 4.56917e-02 * TMath::Log(clusterE) ) / ( 1 + ( 1.10960e+00 * TMath::Exp( ( clusterE - 7.27315e+01 ) / 4.39902e+02 ) ) ) * addCalibFromPi0;
//     else 
//       return 1/2.22*( 1.85096e+00  + 4.56917e-02 * TMath::Log(clusterE) ) / ( 1 + ( 1.10960e+00 * TMath::Exp( ( clusterE - 7.27315e+01 ) / 4.39902e+02 ) ) );
  } else if(caloEnum == kLFHCAL){
//     return 0.9;
    float corrFac1  = (-3.85201e+02+3.86087e+02 * TMath::Erf( (clusterE+4.44436e+01)/(TMath::Sqrt(2)*1.32891e+01)));
    return corrFac1;
//     return tempE;
//     return 1;
  } else if(caloEnum == kHCALOUT){
//     return 0.62;
    float corrFac1  = ( -4.98500e+01+5.05334e+01 * TMath::Erf( (clusterE+6.01498e+01)/(TMath::Sqrt(2)*2.08632e+01)));
    return corrFac1;
    float tempE     = clusterE/corrFac1;
    float corrFac2  = ( 9.74290e-01 - 7.99094e-04 * TMath::Log(tempE) ) / ( 1 + ( -6.64268e-02 * TMath::Exp( ( tempE - 1.32971e+01 ) / -3.23923e+01 ) ) );
    return clusterE /tempE/ corrFac2;
//     return tempE;
  } else if(caloEnum == kHCALIN){
    return 0.62;
  } else if(caloEnum == kEHCAL){
    return 0.37;
  } else {
    return 1;
  }
  return 1.0;
}

float getCorrectionFactorECal(float clusterE, int caloEnum, int algoEnum, int MCtrueID){
  if(caloEnum == kFHCAL){
      return 1.0;
  } else if(caloEnum == kEEMC){
    return 1;
  } else if(caloEnum == kBECAL){
    return 1;
  } else if(caloEnum == kFEMC){
    return 1;
  } else if(caloEnum == kHCALOUT){
    return 1;
  } else if(caloEnum == kHCALIN){
    return 1;
  } else if(caloEnum == kEHCAL){
    return 1;
  } else {
    return 1;
  }
  return 1.0;
}

float getEnergySmearing( int caloEnum, int algoEnum){
  // _fRandom.SetSeed(0);
  // if(caloEnum==kFHCAL){
  //   if(algoEnum==kMA){
  //     return _fRandom.Gaus(1,0.0985);
  //   } else if(algoEnum==kV1){
  //     return _fRandom.Gaus(1,0.10);
  //   } else if(algoEnum==kV3){
  //     return _fRandom.Gaus(1,0.04);
  //   } else {
  //     return 1.0;
  //   }
  // } else {
    return 1.0;
  // }
  // return 1.0;
}

void fillHCalClustersIntoJetFindingInputs( int caloEnum, int clusterizerEnum,  
  std::vector<float> & jetf_hcal_E, std::vector<float> & jetf_hcal_px, std::vector<float> & jetf_hcal_py, std::vector<float> & jetf_hcal_pz,
  std::vector<float> & jetf_calo_E, std::vector<float> & jetf_calo_px, std::vector<float> & jetf_calo_py, std::vector<float> & jetf_calo_pz,
  std::vector<float> & jetf_all_E, std::vector<float> & jetf_all_px, std::vector<float> & jetf_all_py, std::vector<float> & jetf_all_pz
)
{
  for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[clusterizerEnum][caloEnum].size(); iclus++){
      // if(!(_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_isMatched){
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
      // }
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
      // if(!(_clusters_calo[clusterizerEnum][caloEnum].at(iclus)).cluster_isMatched){
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
      // }
  }
}

#endif
