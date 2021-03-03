// ANCHOR basic struct for towers in clusterizer
typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_trueID;
} towersStrct;

enum calotype {
    kFHCAL         = 0,
    kFEMC         = 1,
    kDRCALO        = 2,
    kEEMC         = 3,
    kCEMC         = 4
};

// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }


// conversion functions for processing
float* EtaPhiFromIndices(int ieta,int iphi, int granularity = 1,float energy = 0, int caloSelect = 0);
TVector3 TowerPositionVectorFromIndices(int ieta,int iphi, int caloSelect = 0, int granularity = 1);

float weightM02 = 4.5;
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, Bool_t debugOutput);


// ANCHOR function to return eta and phi from cell iEta and iPhi indices
float * EtaPhiFromIndices(int i_x,int i_y, int granularity = 1, float energy, int caloSelect = 0){
  static float eta_phi[2];
  float zHC = 400;
  float twrD = 10./granularity;
  // center (x=0,y=0) at 27,27 with 10x10 granularity
  int center_x = 27*granularity;
  int center_y = 27*granularity;
  if(caloSelect==kFEMC){
    zHC = 310;
    twrD = 5.535/granularity;
    center_x = 32*granularity;
    center_y = 32*granularity;
  }

  TVector3 vecTwr((center_x-i_x)*twrD,(center_y-i_y)*twrD,zHC);
  eta_phi[0] = vecTwr.Eta();
  eta_phi[1] = (vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());
  if(eta_phi[0]>10)cout << "\t" << center_x-i_x << ":" << center_y-i_y << "\tEta " << eta_phi[0] << "\tPhi "<< eta_phi[1]/**180/TMath::Pi()*/ << "\t" << energy << endl;

  return eta_phi;
}

// ANCHOR function to return a TVector3 for the tower position based on iEta and iPhi indices
TVector3 TowerPositionVectorFromIndices(int i_x,int i_y, int caloSelect = 0, int granularity = 1){
  float zHC = 400;
  float twrD = 10./granularity;
  // center (x=0,y=0) at 27,27 with 10x10 granularity
  int center_x = 27*granularity;
  int center_y = 27*granularity;

  if(caloSelect==kFEMC){
    zHC = 310;
    twrD = 5.535/granularity;
    center_x = 32*granularity;
    center_y = 32*granularity;
  }

  TVector3 twrPositionVec((center_x-i_x)*twrD,(center_y-i_y)*twrD,zHC);
  return twrPositionVec;
}


// ANCHOR function to determine shower shape
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, Bool_t debugOutput){
    static float returnVariables[4]; //0:M02, 1:M20, 2:eta, 3: phi
    float w_tot = 0;
    std::vector<float> w_i;
    TVector3 vecTwr;

    //calculation of weights and weighted position vector
    for(int cellI=0; cellI<cluster_towers.size(); cellI++){
        w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
        w_tot += w_i.at(cellI);
        vecTwr += w_i.at(cellI)*TowerPositionVectorFromIndices(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, caloSelect);
    }
    returnVariables[2]=vecTwr.Eta();
    returnVariables[3]=(vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());

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


// NxN global cluster variables
int _nclusters_NxN_FHCAL = 0;
float* _clusters_NxN_FHCAL_E            = new float[_maxNclusters];
float* _clusters_NxN_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_NxN_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_NxN_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_NxN_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_NxN_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_NxN_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_NxN_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_NxN_FHCAL_NtrueID       = new int[_maxNclusters];

int _nclusters_NxN_FEMC = 0;
float* _clusters_NxN_FEMC_E            = new float[_maxNclusters];
float* _clusters_NxN_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_NxN_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_NxN_FEMC_M02         = new float[_maxNclusters];
float* _clusters_NxN_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_NxN_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_NxN_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_NxN_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_NxN_FEMC_NtrueID       = new int[_maxNclusters];

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

// XN global cluster variables
int _nclusters_XN_FHCAL = 0;
float* _clusters_XN_FHCAL_E            = new float[_maxNclusters];
float* _clusters_XN_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_XN_FHCAL_Phi         = new float[_maxNclusters];
float* _clusters_XN_FHCAL_M02         = new float[_maxNclusters];
float* _clusters_XN_FHCAL_M20         = new float[_maxNclusters];
bool* _clusters_XN_FHCAL_isMatched         = new bool[_maxNclusters];
int* _clusters_XN_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_XN_FHCAL_trueID       = new int[_maxNclusters];
int* _clusters_XN_FHCAL_NtrueID       = new int[_maxNclusters];

int _nclusters_XN_FEMC = 0;
float* _clusters_XN_FEMC_E            = new float[_maxNclusters];
float* _clusters_XN_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_XN_FEMC_Phi         = new float[_maxNclusters];
float* _clusters_XN_FEMC_M02         = new float[_maxNclusters];
float* _clusters_XN_FEMC_M20         = new float[_maxNclusters];
bool* _clusters_XN_FEMC_isMatched         = new bool[_maxNclusters];
int* _clusters_XN_FEMC_NTower         = new int[_maxNclusters];
int* _clusters_XN_FEMC_trueID       = new int[_maxNclusters];
int* _clusters_XN_FEMC_NtrueID       = new int[_maxNclusters];