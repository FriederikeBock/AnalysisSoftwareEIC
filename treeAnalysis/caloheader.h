// ANCHOR basic struct for towers in clusterizer
int _granularityCalo = 1;

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
    kCEMC         = 4
};
TString str_calorimeter[5] = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC"};
const int _active_calo = 2;

enum clusterizertype {
    kV1         = 0,
    kV3         = 1,
    k3x3        = 2,
    k5x5        = 3,
    kC3         = 4,
    kC5         = 5,
    kMA         = 6,
    kDummy      = 7
};
TString str_clusterizer[7] = {"V1", "V3", "3x3", "5x5", "C3", "C5", "MA"};
const int _active_algo = 7;
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

float weightM02 = 4.5;
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, Bool_t debugOutput);


// ANCHOR function to return eta and phi from cell iEta and iPhi indices
float * EtaPhiFromIndices(int i_x,int i_y, float energy, int caloSelect = 0){
  static float eta_phi[2];
  float zHC = 400;
  float twrD = 10./_granularityCalo;
  // center (x=0,y=0) at 27,27 with 10x10 granularity, and 53,53 with 5x5 granularity
  int center_x = _granularityCalo==1 ? 27 : 53;
  int center_y = _granularityCalo==1 ? 27 : 53;
  // center (x=0,y=0) at 27,27 with 5.5x5.5 granularity, and 65,65 with 2.75x2.75 granularity
  if(caloSelect==kFEMC){
    zHC = 310;
    twrD = 5.535/_granularityCalo;
    center_x = _granularityCalo==1 ? 32 : 65;
    center_y = _granularityCalo==1 ? 32 : 65;
  }

  TVector3 vecTwr((center_x-i_x)*twrD,(center_y-i_y)*twrD,zHC);
  eta_phi[0] = vecTwr.Eta();
  eta_phi[1] = (vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());
  if(eta_phi[0]>10)cout << "\t" << center_x-i_x << ":" << center_y-i_y << "\tEta " << eta_phi[0] << "\tPhi "<< eta_phi[1]/**180/TMath::Pi()*/ << "\t" << energy << endl;

  return eta_phi;
}

// ANCHOR function to return a TVector3 for the tower position based on iEta and iPhi indices
TVector3 TowerPositionVectorFromIndices(int i_x,int i_y, int caloSelect = 0){
  float zHC = 350;
  float twrD = 10./_granularityCalo;
  // center (x=0,y=0) at 27,27 with 10x10 granularity
  int center_x = 27*_granularityCalo;
  int center_y = 27*_granularityCalo;

  if(caloSelect==kFEMC){
    zHC = 291;
    twrD = 5.535/_granularityCalo;
    center_x = 32*_granularityCalo;
    center_y = 32*_granularityCalo;
  }

  TVector3 twrPositionVec((i_x-center_x)*twrD,(i_y-center_y)*twrD,zHC);
  return twrPositionVec;
}


// ANCHOR function to determine shower shape
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, bool debugOutput){
    static float returnVariables[7]; //0:M02, 1:M20, 2:eta, 3: phi
    float w_tot = 0;
    std::vector<float> w_i;
    TVector3 vecTwr;
    float zHC = 400;
    if(caloSelect==kFHCAL)zHC=350;
    if(caloSelect==kFEMC)zHC=291;

    vecTwr = {0.,0.,0.};
    //calculation of weights and weighted position vector
    for(int cellI=0; cellI<cluster_towers.size(); cellI++){
        w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
        w_tot += w_i.at(cellI);
        vecTwr += w_i.at(cellI)*TowerPositionVectorFromIndices(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, caloSelect);
    }
    returnVariables[2]=vecTwr.Eta();
    returnVariables[3]=(vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());
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