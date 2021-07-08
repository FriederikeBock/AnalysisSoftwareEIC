// ANCHOR basic struct for towers in clusterizer
typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_trueID;
} towersStrct;

typedef struct {
  int particle_ID;
  float highest_E;
  int nClusters;
} occuranceStrct;

enum calotype {
  kFHCAL         = 0,
  kFEMC         = 1,
  kDRCALO        = 2,
  kEEMC         = 3,
  kCEMC         = 4,
  kEHCAL         = 5,
  kHCALIN       = 6,
  kHCALOUT       = 7,
  kLFHCAL        = 8,
  kEEMCG         = 9,
  kBECAL         = 10
};
const int maxcalo = 11;
int calogeomindex[maxcalo] = {0};

TString str_calorimeter[maxcalo] = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC", "EHCAL", "HCALIN", "HCALOUT", "LFHCAL", "EEMCG", "BECAL"};
const int _active_calo = 11;

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
// TString str_clusterizer[7] = {"V1", "V3", "3x3", "5x5", "C3", "C5", "MA"};
TString str_clusterizer[7] = {"MA", "V3", "V1", "5x5", "C5", "C3", "3x3"};
const int _active_algo = 1;
// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }

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
float* EtaPhiFromIndices(int ieta,int iphi,float energy = 0, int caloSelect = 0);
TVector3 TowerPositionVectorFromIndicesGeometry(int ieta,int iphi, int caloSelect = 0);

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

Int_t ReturnProjectionIndexForCalorimeter(Int_t caloID){
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
    default:
      cout << "ReturnProjectionIndexForCalorimeter: caloID " << caloID << " not defined, returning -1" << endl;
      return -1;
  }
  return -1;
}

// ANCHOR function to return a TVector3 for the tower position based on iEta and iPhi indices
TVector3 TowerPositionVectorFromIndicesGeometry(int i_Eta,int i_Phi, int caloSelect = 0){
  float xpos = -10000;
  float ypos = -10000;
  float zpos = -10000;

  if(calogeomindex[caloSelect]==-1) cout << "calorimeter not found in geometry!" << endl;
  tt_geometry->GetEntry(calogeomindex[caloSelect]);
  for (Long64_t itow=0; itow<_calogeom_towers_N;itow++) {
    if(_calogeom_towers_iEta[itow]==i_Eta){
      if(_calogeom_towers_iPhi[itow]==i_Phi){
        xpos = _calogeom_towers_x[itow];
        ypos = _calogeom_towers_y[itow];
        zpos = _calogeom_towers_z[itow];
      }
    }
  }
  if(xpos==-10000 || ypos==-10000 || zpos==-10000){
    cout << "something is terribly wrong in your geometry: x="<<xpos<<"\ty="<<ypos<<"\tz="<<zpos<<endl;
  }
  TVector3 twrPositionVec(xpos,ypos,zpos);
  return twrPositionVec;
}


// ANCHOR function to determine shower shape
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, bool debugOutput){
    static float returnVariables[7]; //0:M02, 1:M20, 2:eta, 3: phi
    float w_tot = 0;
    std::vector<float> w_i;
    TVector3 vecTwr;
    TVector3 vecTwrTmp;
    float zHC = 400;

    vecTwr = {0.,0.,0.};
    //calculation of weights and weighted position vector
    for(int cellI=0; cellI<cluster_towers.size(); cellI++){
        w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
        w_tot += w_i.at(cellI);
        vecTwrTmp = TowerPositionVectorFromIndicesGeometry(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, caloSelect);
        // cout << caloSelect << "\t" << vecTwrTmp.X() << "\t" << vecTwrTmp.Y() << "\t" << vecTwrTmp.Z() << endl;
        vecTwr += w_i.at(cellI)*vecTwrTmp;
        if(cellI==0)zHC=vecTwrTmp.Z();
    }
    returnVariables[2]=vecTwr.Eta();
    returnVariables[3]=vecTwr.Phi();
    // vecTwr*=1/w_tot;//(zHC/vecTwr.Z());
    vecTwr*=(zHC/vecTwr.Z());
    returnVariables[4]=vecTwr.X();
    returnVariables[5]=vecTwr.Y();
    returnVariables[6]=vecTwr.Z();

    //calculation of M02
    float delta_phi_phi[4] = {0};
    float delta_eta_eta[4] = {0};
    float delta_eta_phi[4] = {0};
    int supModLeadingCell = -1;
    for(int cellI=0; cellI<cluster_towers.size(); cellI++){
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
int _nclusters_3x3_FHCAL = 0;
float* _clusters_3x3_FHCAL_E            = new float[_maxNclusters];
float* _clusters_3x3_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_3x3_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_3x3_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_3x3_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_3x3_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_3x3_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_3x3_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_3x3_FHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_3x3_FHCAL_X         = new float[_maxNclusters];
float* _clusters_3x3_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_3x3_FHCAL_Z         = new float[_maxNclusters];

int _nclusters_3x3_FEMC = 0;
float* _clusters_3x3_FEMC_E            = new float[_maxNclusters];
float* _clusters_3x3_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_3x3_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_3x3_FEMC_M02         = new float[_maxNclusters];
float* _clusters_3x3_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_3x3_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_3x3_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_3x3_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_3x3_FEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_3x3_FEMC_X         = new float[_maxNclusters];
float* _clusters_3x3_FEMC_Y         = new float[_maxNclusters];
float* _clusters_3x3_FEMC_Z         = new float[_maxNclusters];

// 5x5 global cluster variables
int _nclusters_5x5_FHCAL = 0;
float* _clusters_5x5_FHCAL_E            = new float[_maxNclusters];
float* _clusters_5x5_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_5x5_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_5x5_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_5x5_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_5x5_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_5x5_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_5x5_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_5x5_FHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_5x5_FHCAL_X         = new float[_maxNclusters];
float* _clusters_5x5_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_5x5_FHCAL_Z         = new float[_maxNclusters];

int _nclusters_5x5_FEMC = 0;
float* _clusters_5x5_FEMC_E            = new float[_maxNclusters];
float* _clusters_5x5_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_5x5_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_5x5_FEMC_M02         = new float[_maxNclusters];
float* _clusters_5x5_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_5x5_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_5x5_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_5x5_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_5x5_FEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_5x5_FEMC_X         = new float[_maxNclusters];
float* _clusters_5x5_FEMC_Y         = new float[_maxNclusters];
float* _clusters_5x5_FEMC_Z         = new float[_maxNclusters];

// V3 global cluster variables
int _nclusters_V3_FHCAL = 0;
float* _clusters_V3_FHCAL_E            = new float[_maxNclusters];
float* _clusters_V3_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_V3_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_V3_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_V3_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_V3_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_V3_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_V3_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_V3_FHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_V3_FHCAL_X         = new float[_maxNclusters];
float* _clusters_V3_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_V3_FHCAL_Z         = new float[_maxNclusters];

int _nclusters_V3_FEMC = 0;
float* _clusters_V3_FEMC_E            = new float[_maxNclusters];
float* _clusters_V3_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_V3_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_V3_FEMC_M02         = new float[_maxNclusters];
float* _clusters_V3_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_V3_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_V3_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_V3_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_V3_FEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_V3_FEMC_X         = new float[_maxNclusters];
float* _clusters_V3_FEMC_Y         = new float[_maxNclusters];
float* _clusters_V3_FEMC_Z         = new float[_maxNclusters];

// MA global cluster variables
int _nclusters_MA_EEMC = 0;
float* _clusters_MA_EEMC_E            = new float[_maxNclusters];
float* _clusters_MA_EEMC_Eta         = new float[_maxNclusters];
float* _clusters_MA_EEMC_Phi         = new float[_maxNclusters];
float* _clusters_MA_EEMC_M02         = new float[_maxNclusters];
float* _clusters_MA_EEMC_M20         = new float[_maxNclusters];
bool* _clusters_MA_EEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_EEMC_NTower         = new int[_maxNclusters];
int* _clusters_MA_EEMC_trueID       = new int[_maxNclusters];
int* _clusters_MA_EEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_EEMC_X         = new float[_maxNclusters];
float* _clusters_MA_EEMC_Y         = new float[_maxNclusters];
float* _clusters_MA_EEMC_Z         = new float[_maxNclusters];

int _nclusters_MA_EEMCG = 0;
float* _clusters_MA_EEMCG_E            = new float[_maxNclusters];
float* _clusters_MA_EEMCG_Eta         = new float[_maxNclusters];
float* _clusters_MA_EEMCG_Phi         = new float[_maxNclusters];
float* _clusters_MA_EEMCG_M02         = new float[_maxNclusters];
float* _clusters_MA_EEMCG_M20         = new float[_maxNclusters];
bool* _clusters_MA_EEMCG_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_EEMCG_NTower         = new int[_maxNclusters];
int* _clusters_MA_EEMCG_trueID       = new int[_maxNclusters];
int* _clusters_MA_EEMCG_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_EEMCG_X         = new float[_maxNclusters];
float* _clusters_MA_EEMCG_Y         = new float[_maxNclusters];
float* _clusters_MA_EEMCG_Z         = new float[_maxNclusters];

int _nclusters_MA_EHCAL = 0;
float* _clusters_MA_EHCAL_E            = new float[_maxNclusters];
float* _clusters_MA_EHCAL_Eta         = new float[_maxNclusters];
float* _clusters_MA_EHCAL_Phi         = new float[_maxNclusters];
float* _clusters_MA_EHCAL_M02         = new float[_maxNclusters];
float* _clusters_MA_EHCAL_M20         = new float[_maxNclusters];
bool* _clusters_MA_EHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_EHCAL_NTower         = new int[_maxNclusters];
int* _clusters_MA_EHCAL_trueID       = new int[_maxNclusters];
int* _clusters_MA_EHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_EHCAL_X         = new float[_maxNclusters];
float* _clusters_MA_EHCAL_Y         = new float[_maxNclusters];
float* _clusters_MA_EHCAL_Z         = new float[_maxNclusters];

int _nclusters_MA_CEMC = 0;
float* _clusters_MA_CEMC_E            = new float[_maxNclusters];
float* _clusters_MA_CEMC_Eta         = new float[_maxNclusters];
float* _clusters_MA_CEMC_Phi         = new float[_maxNclusters];
float* _clusters_MA_CEMC_M02         = new float[_maxNclusters];
float* _clusters_MA_CEMC_M20         = new float[_maxNclusters];
bool* _clusters_MA_CEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_CEMC_NTower         = new int[_maxNclusters];
int* _clusters_MA_CEMC_trueID       = new int[_maxNclusters];
int* _clusters_MA_CEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_CEMC_X         = new float[_maxNclusters];
float* _clusters_MA_CEMC_Y         = new float[_maxNclusters];
float* _clusters_MA_CEMC_Z         = new float[_maxNclusters];

int _nclusters_MA_BECAL = 0;
float* _clusters_MA_BECAL_E            = new float[_maxNclusters];
float* _clusters_MA_BECAL_Eta         = new float[_maxNclusters];
float* _clusters_MA_BECAL_Phi         = new float[_maxNclusters];
float* _clusters_MA_BECAL_M02         = new float[_maxNclusters];
float* _clusters_MA_BECAL_M20         = new float[_maxNclusters];
bool* _clusters_MA_BECAL_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_BECAL_NTower         = new int[_maxNclusters];
int* _clusters_MA_BECAL_trueID       = new int[_maxNclusters];
int* _clusters_MA_BECAL_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_BECAL_X         = new float[_maxNclusters];
float* _clusters_MA_BECAL_Y         = new float[_maxNclusters];
float* _clusters_MA_BECAL_Z         = new float[_maxNclusters];

int _nclusters_MA_HCALIN = 0;
float* _clusters_MA_HCALIN_E            = new float[_maxNclusters];
float* _clusters_MA_HCALIN_Eta         = new float[_maxNclusters];
float* _clusters_MA_HCALIN_Phi         = new float[_maxNclusters];
float* _clusters_MA_HCALIN_M02         = new float[_maxNclusters];
float* _clusters_MA_HCALIN_M20         = new float[_maxNclusters];
bool* _clusters_MA_HCALIN_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_HCALIN_NTower         = new int[_maxNclusters];
int* _clusters_MA_HCALIN_trueID       = new int[_maxNclusters];
int* _clusters_MA_HCALIN_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_HCALIN_X         = new float[_maxNclusters];
float* _clusters_MA_HCALIN_Y         = new float[_maxNclusters];
float* _clusters_MA_HCALIN_Z         = new float[_maxNclusters];

int _nclusters_MA_HCALOUT = 0;
float* _clusters_MA_HCALOUT_E            = new float[_maxNclusters];
float* _clusters_MA_HCALOUT_Eta         = new float[_maxNclusters];
float* _clusters_MA_HCALOUT_Phi         = new float[_maxNclusters];
float* _clusters_MA_HCALOUT_M02         = new float[_maxNclusters];
float* _clusters_MA_HCALOUT_M20         = new float[_maxNclusters];
bool* _clusters_MA_HCALOUT_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_HCALOUT_NTower         = new int[_maxNclusters];
int* _clusters_MA_HCALOUT_trueID       = new int[_maxNclusters];
int* _clusters_MA_HCALOUT_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_HCALOUT_X         = new float[_maxNclusters];
float* _clusters_MA_HCALOUT_Y         = new float[_maxNclusters];
float* _clusters_MA_HCALOUT_Z         = new float[_maxNclusters];

int _nclusters_MA_LFHCAL = 0;
float* _clusters_MA_LFHCAL_E            = new float[_maxNclusters];
float* _clusters_MA_LFHCAL_Eta         = new float[_maxNclusters];
float* _clusters_MA_LFHCAL_Phi         = new float[_maxNclusters];
float* _clusters_MA_LFHCAL_M02         = new float[_maxNclusters];
float* _clusters_MA_LFHCAL_M20         = new float[_maxNclusters];
bool* _clusters_MA_LFHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_LFHCAL_NTower         = new int[_maxNclusters];
int* _clusters_MA_LFHCAL_trueID       = new int[_maxNclusters];
int* _clusters_MA_LFHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_LFHCAL_X         = new float[_maxNclusters];
float* _clusters_MA_LFHCAL_Y         = new float[_maxNclusters];
float* _clusters_MA_LFHCAL_Z         = new float[_maxNclusters];

int _nclusters_MA_FHCAL = 0;
float* _clusters_MA_FHCAL_E            = new float[_maxNclusters];
float* _clusters_MA_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_MA_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_MA_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_MA_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_MA_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_MA_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_MA_FHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_FHCAL_X         = new float[_maxNclusters];
float* _clusters_MA_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_MA_FHCAL_Z         = new float[_maxNclusters];

int _nclusters_MA_FEMC = 0;
float* _clusters_MA_FEMC_E            = new float[_maxNclusters];
float* _clusters_MA_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_MA_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_MA_FEMC_M02         = new float[_maxNclusters];
float* _clusters_MA_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_MA_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_MA_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_MA_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_MA_FEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_MA_FEMC_X         = new float[_maxNclusters];
float* _clusters_MA_FEMC_Y         = new float[_maxNclusters];
float* _clusters_MA_FEMC_Z         = new float[_maxNclusters];

// C3 global cluster variables
int _nclusters_C3_FHCAL = 0;
float* _clusters_C3_FHCAL_E            = new float[_maxNclusters];
float* _clusters_C3_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_C3_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_C3_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_C3_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_C3_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_C3_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_C3_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_C3_FHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_C3_FHCAL_X         = new float[_maxNclusters];
float* _clusters_C3_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_C3_FHCAL_Z         = new float[_maxNclusters];

int _nclusters_C3_FEMC = 0;
float* _clusters_C3_FEMC_E            = new float[_maxNclusters];
float* _clusters_C3_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_C3_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_C3_FEMC_M02         = new float[_maxNclusters];
float* _clusters_C3_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_C3_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_C3_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_C3_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_C3_FEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_C3_FEMC_X         = new float[_maxNclusters];
float* _clusters_C3_FEMC_Y         = new float[_maxNclusters];
float* _clusters_C3_FEMC_Z         = new float[_maxNclusters];

// C5 global cluster variables
int _nclusters_C5_FHCAL = 0;
float* _clusters_C5_FHCAL_E            = new float[_maxNclusters];
float* _clusters_C5_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_C5_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_C5_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_C5_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_C5_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_C5_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_C5_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_C5_FHCAL_NtrueID       = new int[_maxNclusters];
float* _clusters_C5_FHCAL_X         = new float[_maxNclusters];
float* _clusters_C5_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_C5_FHCAL_Z         = new float[_maxNclusters];

int _nclusters_C5_FEMC = 0;
float* _clusters_C5_FEMC_E            = new float[_maxNclusters];
float* _clusters_C5_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_C5_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_C5_FEMC_M02         = new float[_maxNclusters];
float* _clusters_C5_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_C5_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_C5_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_C5_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_C5_FEMC_NtrueID       = new int[_maxNclusters];
float* _clusters_C5_FEMC_X         = new float[_maxNclusters];
float* _clusters_C5_FEMC_Y         = new float[_maxNclusters];
float* _clusters_C5_FEMC_Z         = new float[_maxNclusters];


// MA global cluster variables
int _nclusters_V1_DRCALO = 0;
float* _clusters_V1_DRCALO_E            = new float[_maxNclusters];
float* _clusters_V1_DRCALO_Eta         = new float[_maxNclusters];
float* _clusters_V1_DRCALO_Phi         = new float[_maxNclusters];
float* _clusters_V1_DRCALO_M02         = new float[_maxNclusters];
float* _clusters_V1_DRCALO_M20         = new float[_maxNclusters];
bool* _clusters_V1_DRCALO_isMatched         = new bool[_maxNclusters];
int* _clusters_V1_DRCALO_NTower         = new int[_maxNclusters];
int* _clusters_V1_DRCALO_trueID       = new int[_maxNclusters];
int* _clusters_V1_DRCALO_NtrueID       = new int[_maxNclusters];
float* _clusters_V1_DRCALO_X         = new float[_maxNclusters];
float* _clusters_V1_DRCALO_Y         = new float[_maxNclusters];
float* _clusters_V1_DRCALO_Z         = new float[_maxNclusters];

bool loadClusterizerInput(
  int clusterizerEnum,
  int caloEnum,
  int &nclusters,
  float* &clusters_E,
  float* &clusters_Eta,
  float* &clusters_Phi,
  float* &clusters_M02,
  float* &clusters_M20,
  bool* &clusters_isMatched,
  int* &clusters_NTower,
  int* &clusters_trueID,
  int* &clusters_NtrueID,
  float* &clusters_X,
  float* &clusters_Y,
  float* &clusters_Z
){
  // cout << clusterizerEnum << "\t" << caloEnum << endl;
  if(clusterizerEnum==k3x3){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_3x3_FHCAL;
      clusters_E = _clusters_3x3_FHCAL_E;
      clusters_Eta = _clusters_3x3_FHCAL_Eta;
      clusters_Phi = _clusters_3x3_FHCAL_Phi;
      clusters_M02 = _clusters_3x3_FHCAL_M02;
      clusters_M20 = _clusters_3x3_FHCAL_M20;
      clusters_isMatched = _clusters_3x3_FHCAL_isMatched;
      clusters_NTower = _clusters_3x3_FHCAL_NTower;
      clusters_trueID = _clusters_3x3_FHCAL_trueID;
      clusters_NtrueID = _clusters_3x3_FHCAL_NtrueID;
      clusters_X = _clusters_3x3_FHCAL_X;
      clusters_Y = _clusters_3x3_FHCAL_Y;
      clusters_Z = _clusters_3x3_FHCAL_Z;
      return true;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_3x3_FEMC;
      clusters_E = _clusters_3x3_FEMC_E;
      clusters_Eta = _clusters_3x3_FEMC_Eta;
      clusters_Phi = _clusters_3x3_FEMC_Phi;
      clusters_M02 = _clusters_3x3_FEMC_M02;
      clusters_M20 = _clusters_3x3_FEMC_M20;
      clusters_isMatched = _clusters_3x3_FEMC_isMatched;
      clusters_NTower = _clusters_3x3_FEMC_NTower;
      clusters_trueID = _clusters_3x3_FEMC_trueID;
      clusters_NtrueID = _clusters_3x3_FEMC_NtrueID;
      clusters_X = _clusters_3x3_FEMC_X;
      clusters_Y = _clusters_3x3_FEMC_Y;
      clusters_Z = _clusters_3x3_FEMC_Z;
      return true;
    } else {
      return false;
      // nclusters = 0;
      // setFLOATClusterArrayToZero(clusters_E);
      // setFLOATClusterArrayToZero(clusters_Eta);
      // setFLOATClusterArrayToZero(clusters_Phi);
      // setFLOATClusterArrayToZero(clusters_M02);
      // setFLOATClusterArrayToZero(clusters_M20);
      // setBOOLClusterArrayToZero(clusters_isMatched);
      // setINTClusterArrayToZero(clusters_NTower);
      // setINTClusterArrayToZero(clusters_trueID);
      // setINTClusterArrayToZero(clusters_NtrueID);
      // setFLOATClusterArrayToZero(clusters_X);
      // setFLOATClusterArrayToZero(clusters_Y);
      // setFLOATClusterArrayToZero(clusters_Z);
    }
  }
  else if(clusterizerEnum==k5x5){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_5x5_FHCAL;
      clusters_E = _clusters_5x5_FHCAL_E;
      clusters_Eta = _clusters_5x5_FHCAL_Eta;
      clusters_Phi = _clusters_5x5_FHCAL_Phi;
      clusters_M02 = _clusters_5x5_FHCAL_M02;
      clusters_M20 = _clusters_5x5_FHCAL_M20;
      clusters_isMatched = _clusters_5x5_FHCAL_isMatched;
      clusters_NTower = _clusters_5x5_FHCAL_NTower;
      clusters_trueID = _clusters_5x5_FHCAL_trueID;
      clusters_NtrueID = _clusters_5x5_FHCAL_NtrueID;
      clusters_X = _clusters_5x5_FHCAL_X;
      clusters_Y = _clusters_5x5_FHCAL_Y;
      clusters_Z = _clusters_5x5_FHCAL_Z;
      return true;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_5x5_FEMC;
      clusters_E = _clusters_5x5_FEMC_E;
      clusters_Eta = _clusters_5x5_FEMC_Eta;
      clusters_Phi = _clusters_5x5_FEMC_Phi;
      clusters_M02 = _clusters_5x5_FEMC_M02;
      clusters_M20 = _clusters_5x5_FEMC_M20;
      clusters_isMatched = _clusters_5x5_FEMC_isMatched;
      clusters_NTower = _clusters_5x5_FEMC_NTower;
      clusters_trueID = _clusters_5x5_FEMC_trueID;
      clusters_NtrueID = _clusters_5x5_FEMC_NtrueID;
      clusters_X = _clusters_5x5_FEMC_X;
      clusters_Y = _clusters_5x5_FEMC_Y;
      clusters_Z = _clusters_5x5_FEMC_Z;
      return true;
    } else {
      return false;
      // nclusters = 0;
      // setFLOATClusterArrayToZero(clusters_E);
      // setFLOATClusterArrayToZero(clusters_Eta);
      // setFLOATClusterArrayToZero(clusters_Phi);
      // setFLOATClusterArrayToZero(clusters_M02);
      // setFLOATClusterArrayToZero(clusters_M20);
      // setBOOLClusterArrayToZero(clusters_isMatched);
      // setINTClusterArrayToZero(clusters_NTower);
      // setINTClusterArrayToZero(clusters_trueID);
      // setINTClusterArrayToZero(clusters_NtrueID);
      // setFLOATClusterArrayToZero(clusters_X);
      // setFLOATClusterArrayToZero(clusters_Y);
      // setFLOATClusterArrayToZero(clusters_Z);
    }
  }
// V3 global cluster variables
  else if(clusterizerEnum==kV3){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_V3_FHCAL;
      clusters_E = _clusters_V3_FHCAL_E;
      clusters_Eta = _clusters_V3_FHCAL_Eta;
      clusters_Phi = _clusters_V3_FHCAL_Phi;
      clusters_M02 = _clusters_V3_FHCAL_M02;
      clusters_M20 = _clusters_V3_FHCAL_M20;
      clusters_isMatched = _clusters_V3_FHCAL_isMatched;
      clusters_NTower = _clusters_V3_FHCAL_NTower;
      clusters_trueID = _clusters_V3_FHCAL_trueID;
      clusters_NtrueID = _clusters_V3_FHCAL_NtrueID;
      clusters_X = _clusters_V3_FHCAL_X;
      clusters_Y = _clusters_V3_FHCAL_Y;
      clusters_Z = _clusters_V3_FHCAL_Z;
      return true;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_V3_FEMC;
      clusters_E = _clusters_V3_FEMC_E;
      clusters_Eta = _clusters_V3_FEMC_Eta;
      clusters_Phi = _clusters_V3_FEMC_Phi;
      clusters_M02 = _clusters_V3_FEMC_M02;
      clusters_M20 = _clusters_V3_FEMC_M20;
      clusters_isMatched = _clusters_V3_FEMC_isMatched;
      clusters_NTower = _clusters_V3_FEMC_NTower;
      clusters_trueID = _clusters_V3_FEMC_trueID;
      clusters_NtrueID = _clusters_V3_FEMC_NtrueID;
      clusters_X = _clusters_V3_FEMC_X;
      clusters_Y = _clusters_V3_FEMC_Y;
      clusters_Z = _clusters_V3_FEMC_Z;
      return true;
    } else {
      return false;
      // nclusters = 0;
      // setFLOATClusterArrayToZero(clusters_E);
      // setFLOATClusterArrayToZero(clusters_Eta);
      // setFLOATClusterArrayToZero(clusters_Phi);
      // setFLOATClusterArrayToZero(clusters_M02);
      // setFLOATClusterArrayToZero(clusters_M20);
      // setBOOLClusterArrayToZero(clusters_isMatched);
      // setINTClusterArrayToZero(clusters_NTower);
      // setINTClusterArrayToZero(clusters_trueID);
      // setINTClusterArrayToZero(clusters_NtrueID);
      // setFLOATClusterArrayToZero(clusters_X);
      // setFLOATClusterArrayToZero(clusters_Y);
      // setFLOATClusterArrayToZero(clusters_Z);
    }
  }
// MA global cluster variables
  else if(clusterizerEnum==kMA){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_MA_FHCAL;
      clusters_E = _clusters_MA_FHCAL_E;
      clusters_Eta = _clusters_MA_FHCAL_Eta;
      clusters_Phi = _clusters_MA_FHCAL_Phi;
      clusters_M02 = _clusters_MA_FHCAL_M02;
      clusters_M20 = _clusters_MA_FHCAL_M20;
      clusters_isMatched = _clusters_MA_FHCAL_isMatched;
      clusters_NTower = _clusters_MA_FHCAL_NTower;
      clusters_trueID = _clusters_MA_FHCAL_trueID;
      clusters_NtrueID = _clusters_MA_FHCAL_NtrueID;
      clusters_X = _clusters_MA_FHCAL_X;
      clusters_Y = _clusters_MA_FHCAL_Y;
      clusters_Z = _clusters_MA_FHCAL_Z;
      return true;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_MA_FEMC;
      clusters_E = _clusters_MA_FEMC_E;
      clusters_Eta = _clusters_MA_FEMC_Eta;
      clusters_Phi = _clusters_MA_FEMC_Phi;
      clusters_M02 = _clusters_MA_FEMC_M02;
      clusters_M20 = _clusters_MA_FEMC_M20;
      clusters_isMatched = _clusters_MA_FEMC_isMatched;
      clusters_NTower = _clusters_MA_FEMC_NTower;
      clusters_trueID = _clusters_MA_FEMC_trueID;
      clusters_NtrueID = _clusters_MA_FEMC_NtrueID;
      clusters_X = _clusters_MA_FEMC_X;
      clusters_Y = _clusters_MA_FEMC_Y;
      clusters_Z = _clusters_MA_FEMC_Z;
      return true;
    } else if(caloEnum==kLFHCAL){
      nclusters = _nclusters_MA_LFHCAL;
      clusters_E = _clusters_MA_LFHCAL_E;
      clusters_Eta = _clusters_MA_LFHCAL_Eta;
      clusters_Phi = _clusters_MA_LFHCAL_Phi;
      clusters_M02 = _clusters_MA_LFHCAL_M02;
      clusters_M20 = _clusters_MA_LFHCAL_M20;
      clusters_isMatched = _clusters_MA_LFHCAL_isMatched;
      clusters_NTower = _clusters_MA_LFHCAL_NTower;
      clusters_trueID = _clusters_MA_LFHCAL_trueID;
      clusters_NtrueID = _clusters_MA_LFHCAL_NtrueID;
      clusters_X = _clusters_MA_LFHCAL_X;
      clusters_Y = _clusters_MA_LFHCAL_Y;
      clusters_Z = _clusters_MA_LFHCAL_Z;
      return true;
    } else if(caloEnum==kEEMC){
      nclusters = _nclusters_MA_EEMC;
      clusters_E = _clusters_MA_EEMC_E;
      clusters_Eta = _clusters_MA_EEMC_Eta;
      clusters_Phi = _clusters_MA_EEMC_Phi;
      clusters_M02 = _clusters_MA_EEMC_M02;
      clusters_M20 = _clusters_MA_EEMC_M20;
      clusters_isMatched = _clusters_MA_EEMC_isMatched;
      clusters_NTower = _clusters_MA_EEMC_NTower;
      clusters_trueID = _clusters_MA_EEMC_trueID;
      clusters_NtrueID = _clusters_MA_EEMC_NtrueID;
      clusters_X = _clusters_MA_EEMC_X;
      clusters_Y = _clusters_MA_EEMC_Y;
      clusters_Z = _clusters_MA_EEMC_Z;
      return true;
    } else if(caloEnum==kEEMCG){
      nclusters = _nclusters_MA_EEMCG;
      clusters_E = _clusters_MA_EEMCG_E;
      clusters_Eta = _clusters_MA_EEMCG_Eta;
      clusters_Phi = _clusters_MA_EEMCG_Phi;
      clusters_M02 = _clusters_MA_EEMCG_M02;
      clusters_M20 = _clusters_MA_EEMCG_M20;
      clusters_isMatched = _clusters_MA_EEMCG_isMatched;
      clusters_NTower = _clusters_MA_EEMCG_NTower;
      clusters_trueID = _clusters_MA_EEMCG_trueID;
      clusters_NtrueID = _clusters_MA_EEMCG_NtrueID;
      clusters_X = _clusters_MA_EEMCG_X;
      clusters_Y = _clusters_MA_EEMCG_Y;
      clusters_Z = _clusters_MA_EEMCG_Z;
      return true;
    } else if(caloEnum==kEHCAL){
      nclusters = _nclusters_MA_EHCAL;
      clusters_E = _clusters_MA_EHCAL_E;
      clusters_Eta = _clusters_MA_EHCAL_Eta;
      clusters_Phi = _clusters_MA_EHCAL_Phi;
      clusters_M02 = _clusters_MA_EHCAL_M02;
      clusters_M20 = _clusters_MA_EHCAL_M20;
      clusters_isMatched = _clusters_MA_EHCAL_isMatched;
      clusters_NTower = _clusters_MA_EHCAL_NTower;
      clusters_trueID = _clusters_MA_EHCAL_trueID;
      clusters_NtrueID = _clusters_MA_EHCAL_NtrueID;
      clusters_X = _clusters_MA_EHCAL_X;
      clusters_Y = _clusters_MA_EHCAL_Y;
      clusters_Z = _clusters_MA_EHCAL_Z;
      return true;
    } else if(caloEnum==kCEMC){
      nclusters = _nclusters_MA_CEMC;
      clusters_E = _clusters_MA_CEMC_E;
      clusters_Eta = _clusters_MA_CEMC_Eta;
      clusters_Phi = _clusters_MA_CEMC_Phi;
      clusters_M02 = _clusters_MA_CEMC_M02;
      clusters_M20 = _clusters_MA_CEMC_M20;
      clusters_isMatched = _clusters_MA_CEMC_isMatched;
      clusters_NTower = _clusters_MA_CEMC_NTower;
      clusters_trueID = _clusters_MA_CEMC_trueID;
      clusters_NtrueID = _clusters_MA_CEMC_NtrueID;
      clusters_X = _clusters_MA_CEMC_X;
      clusters_Y = _clusters_MA_CEMC_Y;
      clusters_Z = _clusters_MA_CEMC_Z;
      return true;
    } else if(caloEnum==kHCALIN){
      nclusters = _nclusters_MA_HCALIN;
      clusters_E = _clusters_MA_HCALIN_E;
      clusters_Eta = _clusters_MA_HCALIN_Eta;
      clusters_Phi = _clusters_MA_HCALIN_Phi;
      clusters_M02 = _clusters_MA_HCALIN_M02;
      clusters_M20 = _clusters_MA_HCALIN_M20;
      clusters_isMatched = _clusters_MA_HCALIN_isMatched;
      clusters_NTower = _clusters_MA_HCALIN_NTower;
      clusters_trueID = _clusters_MA_HCALIN_trueID;
      clusters_NtrueID = _clusters_MA_HCALIN_NtrueID;
      clusters_X = _clusters_MA_HCALIN_X;
      clusters_Y = _clusters_MA_HCALIN_Y;
      clusters_Z = _clusters_MA_HCALIN_Z;
      return true;
    } else if(caloEnum==kHCALOUT){
      nclusters = _nclusters_MA_HCALOUT;
      clusters_E = _clusters_MA_HCALOUT_E;
      clusters_Eta = _clusters_MA_HCALOUT_Eta;
      clusters_Phi = _clusters_MA_HCALOUT_Phi;
      clusters_M02 = _clusters_MA_HCALOUT_M02;
      clusters_M20 = _clusters_MA_HCALOUT_M20;
      clusters_isMatched = _clusters_MA_HCALOUT_isMatched;
      clusters_NTower = _clusters_MA_HCALOUT_NTower;
      clusters_trueID = _clusters_MA_HCALOUT_trueID;
      clusters_NtrueID = _clusters_MA_HCALOUT_NtrueID;
      clusters_X = _clusters_MA_HCALOUT_X;
      clusters_Y = _clusters_MA_HCALOUT_Y;
      clusters_Z = _clusters_MA_HCALOUT_Z;
      return true;
    } else {
      return false;
      // nclusters = 0;
      // setFLOATClusterArrayToZero(clusters_E);
      // setFLOATClusterArrayToZero(clusters_Eta);
      // setFLOATClusterArrayToZero(clusters_Phi);
      // setFLOATClusterArrayToZero(clusters_M02);
      // setFLOATClusterArrayToZero(clusters_M20);
      // setBOOLClusterArrayToZero(clusters_isMatched);
      // setINTClusterArrayToZero(clusters_NTower);
      // setINTClusterArrayToZero(clusters_trueID);
      // setINTClusterArrayToZero(clusters_NtrueID);
      // setFLOATClusterArrayToZero(clusters_X);
      // setFLOATClusterArrayToZero(clusters_Y);
      // setFLOATClusterArrayToZero(clusters_Z);
    }
  }
// V1 global cluster variables
  else if(clusterizerEnum==kV1){
    if(caloEnum==kDRCALO){
      nclusters = _nclusters_V1_DRCALO;
      clusters_E = _clusters_V1_DRCALO_E;
      clusters_Eta = _clusters_V1_DRCALO_Eta;
      clusters_Phi = _clusters_V1_DRCALO_Phi;
      clusters_M02 = _clusters_V1_DRCALO_M02;
      clusters_M20 = _clusters_V1_DRCALO_M20;
      clusters_isMatched = _clusters_V1_DRCALO_isMatched;
      clusters_NTower = _clusters_V1_DRCALO_NTower;
      clusters_trueID = _clusters_V1_DRCALO_trueID;
      clusters_NtrueID = _clusters_V1_DRCALO_NtrueID;
      clusters_X = _clusters_V1_DRCALO_X;
      clusters_Y = _clusters_V1_DRCALO_Y;
      clusters_Z = _clusters_V1_DRCALO_Z;
      return true;
    }else if(caloEnum==kFHCAL){
      nclusters = _nclusters_FHCAL;
      clusters_E = _clusters_FHCAL_E;
      clusters_Eta = _clusters_FHCAL_Eta;
      clusters_Phi = _clusters_FHCAL_Phi;
      clusters_M02 = new float[_maxNclusters];
      clusters_M20 = new float[_maxNclusters];
      clusters_isMatched = new bool[_maxNclusters];
      setFLOATClusterArrayToZero(clusters_M02);
      setFLOATClusterArrayToZero(clusters_M20);
      setBOOLClusterArrayToZero(clusters_isMatched);
      clusters_NTower = _clusters_FHCAL_NTower;
      clusters_trueID = _clusters_FHCAL_trueID;
      clusters_NtrueID = new int[_maxNclusters];
      clusters_X = new float[_maxNclusters];
      clusters_Y = new float[_maxNclusters];
      clusters_Z = new float[_maxNclusters];
      setINTClusterArrayToZero(clusters_NtrueID);
      setFLOATClusterArrayToZero(clusters_X);
      setFLOATClusterArrayToZero(clusters_Y);
      setFLOATClusterArrayToZero(clusters_Z);
      return true;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_FEMC;
      clusters_E = _clusters_FEMC_E;
      clusters_Eta = _clusters_FEMC_Eta;
      clusters_Phi = _clusters_FEMC_Phi;
      clusters_M02 = new float[_maxNclusters];
      clusters_M20 = new float[_maxNclusters];
      clusters_isMatched = new bool[_maxNclusters];
      setFLOATClusterArrayToZero(clusters_M02);
      setFLOATClusterArrayToZero(clusters_M20);
      setBOOLClusterArrayToZero(clusters_isMatched);
      clusters_NTower = _clusters_FEMC_NTower;
      clusters_trueID = _clusters_FEMC_trueID;
      clusters_NtrueID = new int[_maxNclusters];
      clusters_X = new float[_maxNclusters];
      clusters_Y = new float[_maxNclusters];
      clusters_Z = new float[_maxNclusters];
      setINTClusterArrayToZero(clusters_NtrueID);
      setFLOATClusterArrayToZero(clusters_X);
      setFLOATClusterArrayToZero(clusters_Y);
      setFLOATClusterArrayToZero(clusters_Z);
      return true;
    } else {
      return false;
      // nclusters = 0;
      // setFLOATClusterArrayToZero(clusters_E);
      // setFLOATClusterArrayToZero(clusters_Eta);
      // setFLOATClusterArrayToZero(clusters_Phi);
      // setFLOATClusterArrayToZero(clusters_M02);
      // setFLOATClusterArrayToZero(clusters_M20);
      // setBOOLClusterArrayToZero(clusters_isMatched);
      // setINTClusterArrayToZero(clusters_NTower);
      // setINTClusterArrayToZero(clusters_trueID);
      // setINTClusterArrayToZero(clusters_NtrueID);
      // setFLOATClusterArrayToZero(clusters_X);
      // setFLOATClusterArrayToZero(clusters_Y);
      // setFLOATClusterArrayToZero(clusters_Z);
    }
  }

// C3 global cluster variables
  else if(clusterizerEnum==kC3){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_C3_FHCAL;
      clusters_E = _clusters_C3_FHCAL_E;
      clusters_Eta = _clusters_C3_FHCAL_Eta;
      clusters_Phi = _clusters_C3_FHCAL_Phi;
      clusters_M02 = _clusters_C3_FHCAL_M02;
      clusters_M20 = _clusters_C3_FHCAL_M20;
      clusters_isMatched = _clusters_C3_FHCAL_isMatched;
      clusters_NTower = _clusters_C3_FHCAL_NTower;
      clusters_trueID = _clusters_C3_FHCAL_trueID;
      clusters_NtrueID = _clusters_C3_FHCAL_NtrueID;
      clusters_X = _clusters_C3_FHCAL_X;
      clusters_Y = _clusters_C3_FHCAL_Y;
      clusters_Z = _clusters_C3_FHCAL_Z;
      return true;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_C3_FEMC;
      clusters_E = _clusters_C3_FEMC_E;
      clusters_Eta = _clusters_C3_FEMC_Eta;
      clusters_Phi = _clusters_C3_FEMC_Phi;
      clusters_M02 = _clusters_C3_FEMC_M02;
      clusters_M20 = _clusters_C3_FEMC_M20;
      clusters_isMatched = _clusters_C3_FEMC_isMatched;
      clusters_NTower = _clusters_C3_FEMC_NTower;
      clusters_trueID = _clusters_C3_FEMC_trueID;
      clusters_NtrueID = _clusters_C3_FEMC_NtrueID;
      clusters_X = _clusters_C3_FEMC_X;
      clusters_Y = _clusters_C3_FEMC_Y;
      clusters_Z = _clusters_C3_FEMC_Z;
      return true;
    } else {
      return false;
      // nclusters = 0;
      // setFLOATClusterArrayToZero(clusters_E);
      // setFLOATClusterArrayToZero(clusters_Eta);
      // setFLOATClusterArrayToZero(clusters_Phi);
      // setFLOATClusterArrayToZero(clusters_M02);
      // setFLOATClusterArrayToZero(clusters_M20);
      // setBOOLClusterArrayToZero(clusters_isMatched);
      // setINTClusterArrayToZero(clusters_NTower);
      // setINTClusterArrayToZero(clusters_trueID);
      // setINTClusterArrayToZero(clusters_NtrueID);
      // setFLOATClusterArrayToZero(clusters_X);
      // setFLOATClusterArrayToZero(clusters_Y);
      // setFLOATClusterArrayToZero(clusters_Z);
    }
  }

// C5 global cluster variables
  else if(clusterizerEnum==kC5){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_C5_FHCAL;
      clusters_E = _clusters_C5_FHCAL_E;
      clusters_Eta = _clusters_C5_FHCAL_Eta;
      clusters_Phi = _clusters_C5_FHCAL_Phi;
      clusters_M02 = _clusters_C5_FHCAL_M02;
      clusters_M20 = _clusters_C5_FHCAL_M20;
      clusters_isMatched = _clusters_C5_FHCAL_isMatched;
      clusters_NTower = _clusters_C5_FHCAL_NTower;
      clusters_trueID = _clusters_C5_FHCAL_trueID;
      clusters_NtrueID = _clusters_C5_FHCAL_NtrueID;
      clusters_X = _clusters_C5_FHCAL_X;
      clusters_Y = _clusters_C5_FHCAL_Y;
      clusters_Z = _clusters_C5_FHCAL_Z;
      return true;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_C5_FEMC;
      clusters_E = _clusters_C5_FEMC_E;
      clusters_Eta = _clusters_C5_FEMC_Eta;
      clusters_Phi = _clusters_C5_FEMC_Phi;
      clusters_M02 = _clusters_C5_FEMC_M02;
      clusters_M20 = _clusters_C5_FEMC_M20;
      clusters_isMatched = _clusters_C5_FEMC_isMatched;
      clusters_NTower = _clusters_C5_FEMC_NTower;
      clusters_trueID = _clusters_C5_FEMC_trueID;
      clusters_NtrueID = _clusters_C5_FEMC_NtrueID;
      clusters_X = _clusters_C5_FEMC_X;
      clusters_Y = _clusters_C5_FEMC_Y;
      clusters_Z = _clusters_C5_FEMC_Z;
      return true;
    } else {
      return false;
      // nclusters = 0;
      // setFLOATClusterArrayToZero(clusters_E);
      // setFLOATClusterArrayToZero(clusters_Eta);
      // setFLOATClusterArrayToZero(clusters_Phi);
      // setFLOATClusterArrayToZero(clusters_M02);
      // setFLOATClusterArrayToZero(clusters_M20);
      // setBOOLClusterArrayToZero(clusters_isMatched);
      // setINTClusterArrayToZero(clusters_NTower);
      // setINTClusterArrayToZero(clusters_trueID);
      // setINTClusterArrayToZero(clusters_NtrueID);
      // setFLOATClusterArrayToZero(clusters_X);
      // setFLOATClusterArrayToZero(clusters_Y);
      // setFLOATClusterArrayToZero(clusters_Z);
    }
  }
  else {
    return false;
    // nclusters = 0;
    // setFLOATClusterArrayToZero(clusters_E);
    // setFLOATClusterArrayToZero(clusters_Eta);
    // setFLOATClusterArrayToZero(clusters_Phi);
    // setFLOATClusterArrayToZero(clusters_M02);
    // setFLOATClusterArrayToZero(clusters_M20);
    // setBOOLClusterArrayToZero(clusters_isMatched);
    // setINTClusterArrayToZero(clusters_NTower);
    // setINTClusterArrayToZero(clusters_trueID);
    // setINTClusterArrayToZero(clusters_NtrueID);
    // setFLOATClusterArrayToZero(clusters_X);
    // setFLOATClusterArrayToZero(clusters_Y);
    // setFLOATClusterArrayToZero(clusters_Z);
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
