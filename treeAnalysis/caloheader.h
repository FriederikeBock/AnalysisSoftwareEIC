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
TString str_calorimeter[5] = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC"};

enum clusterizertype {
    kV1         = 0,
    kV3         = 1,
    kNxN        = 2,
    kXN         = 3
};
TString str_clusterizer[4] = {"V1", "V3", "NxN", "XN"};

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
  float zHC = 350;
  float twrD = 10./granularity;
  // center (x=0,y=0) at 27,27 with 10x10 granularity
  int center_x = 27*granularity;
  int center_y = 27*granularity;

  if(caloSelect==kFEMC){
    zHC = 291;
    twrD = 5.535/granularity;
    center_x = 32*granularity;
    center_y = 32*granularity;
  }

  TVector3 twrPositionVec((center_x-i_x)*twrD,(center_y-i_y)*twrD,zHC);
  return twrPositionVec;
}


// ANCHOR function to determine shower shape
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, int caloSelect, Bool_t debugOutput){
    static float returnVariables[7]; //0:M02, 1:M20, 2:eta, 3: phi
    float w_tot = 0;
    std::vector<float> w_i;
    TVector3 vecTwr;
    float zHC = 400;
    if(caloSelect==kFHCAL)zHC=350;
    if(caloSelect==kFEMC)zHC=291;

    //calculation of weights and weighted position vector
    for(int cellI=0; cellI<cluster_towers.size(); cellI++){
        w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
        w_tot += w_i.at(cellI);
        vecTwr += w_i.at(cellI)*TowerPositionVectorFromIndices(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, caloSelect);
    }
    returnVariables[2]=vecTwr.Eta();
    returnVariables[3]=(vecTwr.Phi()<0 ? vecTwr.Phi()+TMath::Pi() : vecTwr.Phi()-TMath::Pi());
    vecTwr*=(zHC/vecTwr.Z());
    // returnVariables[4]=TowerPositionVectorFromIndices(cluster_towers.at(0).tower_iEta,cluster_towers.at(0).tower_iPhi, caloSelect).X();
    // returnVariables[5]=TowerPositionVectorFromIndices(cluster_towers.at(0).tower_iEta,cluster_towers.at(0).tower_iPhi, caloSelect).Y();
    // returnVariables[6]=TowerPositionVectorFromIndices(cluster_towers.at(0).tower_iEta,cluster_towers.at(0).tower_iPhi, caloSelect).Z();
    returnVariables[4]=vecTwr.X();
    returnVariables[5]=vecTwr.Y();
    returnVariables[6]=vecTwr.Z();
    // if(vecTwr.X()>1e3 || (abs(vecTwr.X()) < 1 && abs(vecTwr.X()) !=0)) cout << "CALO " << caloSelect << "\tx: " << vecTwr.X()<< "\ty: " << vecTwr.Y()<< "\tz: " << vecTwr.Z() << endl;
    if((vecTwr.Z()>351) || (vecTwr.Z()<290)){
      cout << "CALO " << caloSelect << "\tx: " << vecTwr.X()<< "\ty: " << vecTwr.Y()<< "\tz: " << vecTwr.Z() << endl;
      for(int cellI=0; cellI<cluster_towers.size(); cellI++){
        TVector3 vectemp = TowerPositionVectorFromIndices(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, caloSelect);
        cout << "\tx: " << vectemp.X()<< "\ty: " << vectemp.Y()<< "\tz: " << vectemp.Z() << endl;
      }
    }

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
float* _clusters_NxN_FHCAL_X         = new float[_maxNclusters];
float* _clusters_NxN_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_NxN_FHCAL_Z         = new float[_maxNclusters];

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
float* _clusters_NxN_FEMC_X         = new float[_maxNclusters];
float* _clusters_NxN_FEMC_Y         = new float[_maxNclusters];
float* _clusters_NxN_FEMC_Z         = new float[_maxNclusters];

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
float* _clusters_XN_FHCAL_X         = new float[_maxNclusters];
float* _clusters_XN_FHCAL_Y         = new float[_maxNclusters];
float* _clusters_XN_FHCAL_Z         = new float[_maxNclusters];

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
float* _clusters_XN_FEMC_X         = new float[_maxNclusters];
float* _clusters_XN_FEMC_Y         = new float[_maxNclusters];
float* _clusters_XN_FEMC_Z         = new float[_maxNclusters];


void loadClusterizerInput(
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
  if(clusterizerEnum==kNxN){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_NxN_FHCAL;
      clusters_E = _clusters_NxN_FHCAL_E;
      clusters_Eta = _clusters_NxN_FHCAL_Eta;
      clusters_Phi = _clusters_NxN_FHCAL_Phi;
      clusters_M02 = _clusters_NxN_FHCAL_M02;
      clusters_M20 = _clusters_NxN_FHCAL_M20;
      clusters_isMatched = _clusters_NxN_FHCAL_isMatched;
      clusters_NTower = _clusters_NxN_FHCAL_NTower;
      clusters_trueID = _clusters_NxN_FHCAL_trueID;
      clusters_NtrueID = _clusters_NxN_FHCAL_NtrueID;
      clusters_X = _clusters_NxN_FHCAL_X;
      clusters_Y = _clusters_NxN_FHCAL_Y;
      clusters_Z = _clusters_NxN_FHCAL_Z;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_NxN_FEMC;
      clusters_E = _clusters_NxN_FEMC_E;
      clusters_Eta = _clusters_NxN_FEMC_Eta;
      clusters_Phi = _clusters_NxN_FEMC_Phi;
      clusters_M02 = _clusters_NxN_FEMC_M02;
      clusters_M20 = _clusters_NxN_FEMC_M20;
      clusters_isMatched = _clusters_NxN_FEMC_isMatched;
      clusters_NTower = _clusters_NxN_FEMC_NTower;
      clusters_trueID = _clusters_NxN_FEMC_trueID;
      clusters_NtrueID = _clusters_NxN_FEMC_NtrueID;
      clusters_X = _clusters_NxN_FEMC_X;
      clusters_Y = _clusters_NxN_FEMC_Y;
      clusters_Z = _clusters_NxN_FEMC_Z;
    } else {
      nclusters = 0;
      clusters_E = {0};
      clusters_Eta = {0};
      clusters_Phi = {0};
      clusters_M02 = {0};
      clusters_M20 = {0};
      clusters_isMatched = {0};
      clusters_NTower = {0};
      clusters_trueID = {0};
      clusters_NtrueID = {0};
      clusters_X = {0};
      clusters_Y = {0};
      clusters_Z = {0};
    }
  }
// V3 global cluster variables
  if(clusterizerEnum==kV3){
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
    } else {
      nclusters = 0;
      clusters_E = {0};
      clusters_Eta = {0};
      clusters_Phi = {0};
      clusters_M02 = {0};
      clusters_M20 = {0};
      clusters_isMatched = {0};
      clusters_NTower = {0};
      clusters_trueID = {0};
      clusters_NtrueID = {0};
      clusters_X = {0};
      clusters_Y = {0};
      clusters_Z = {0};
    }
  }
// V1 global cluster variables
  if(clusterizerEnum==kV1){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_FHCAL;
      clusters_E = _clusters_FHCAL_E;
      clusters_Eta = _clusters_FHCAL_Eta;
      clusters_Phi = _clusters_FHCAL_Phi;
      clusters_M02 = {0};
      clusters_M20 = {0};
      clusters_isMatched = {0};
      clusters_NTower = _clusters_FHCAL_NTower;
      clusters_trueID = _clusters_FHCAL_trueID;
      clusters_NtrueID = {0};
      clusters_X = {0};
      clusters_Y = {0};
      clusters_Z = {0};
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_FEMC;
      clusters_E = _clusters_FEMC_E;
      clusters_Eta = _clusters_FEMC_Eta;
      clusters_Phi = _clusters_FEMC_Phi;
      clusters_M02 = {0};
      clusters_M20 = {0};
      clusters_isMatched = {0};
      clusters_NTower = _clusters_FEMC_NTower;
      clusters_trueID = _clusters_FEMC_trueID;
      clusters_NtrueID = {0};
      clusters_X = {0};
      clusters_Y = {0};
      clusters_Z = {0};
    } else {
      nclusters = 0;
      clusters_E = {0};
      clusters_Eta = {0};
      clusters_Phi = {0};
      clusters_M02 = {0};
      clusters_M20 = {0};
      clusters_isMatched = {0};
      clusters_NTower = {0};
      clusters_trueID = {0};
      clusters_NtrueID = {0};
      clusters_X = {0};
      clusters_Y = {0};
      clusters_Z = {0};
    }
  }

// XN global cluster variables
  if(clusterizerEnum==kXN){
    if(caloEnum==kFHCAL){
      nclusters = _nclusters_XN_FHCAL;
      clusters_E = _clusters_XN_FHCAL_E;
      clusters_Eta = _clusters_XN_FHCAL_Eta;
      clusters_Phi = _clusters_XN_FHCAL_Phi;
      clusters_M02 = _clusters_XN_FHCAL_M02;
      clusters_M20 = _clusters_XN_FHCAL_M20;
      clusters_isMatched = _clusters_XN_FHCAL_isMatched;
      clusters_NTower = _clusters_XN_FHCAL_NTower;
      clusters_trueID = _clusters_XN_FHCAL_trueID;
      clusters_NtrueID = _clusters_XN_FHCAL_NtrueID;
      clusters_X = _clusters_XN_FHCAL_X;
      clusters_Y = _clusters_XN_FHCAL_Y;
      clusters_Z = _clusters_XN_FHCAL_Z;
    } else if(caloEnum==kFEMC){
      nclusters = _nclusters_XN_FEMC;
      clusters_E = _clusters_XN_FEMC_E;
      clusters_Eta = _clusters_XN_FEMC_Eta;
      clusters_Phi = _clusters_XN_FEMC_Phi;
      clusters_M02 = _clusters_XN_FEMC_M02;
      clusters_M20 = _clusters_XN_FEMC_M20;
      clusters_isMatched = _clusters_XN_FEMC_isMatched;
      clusters_NTower = _clusters_XN_FEMC_NTower;
      clusters_trueID = _clusters_XN_FEMC_trueID;
      clusters_NtrueID = _clusters_XN_FEMC_NtrueID;
      clusters_X = _clusters_XN_FEMC_X;
      clusters_Y = _clusters_XN_FEMC_Y;
      clusters_Z = _clusters_XN_FEMC_Z;
    } else {
      nclusters = 0;
      clusters_E = {0};
      clusters_Eta = {0};
      clusters_Phi = {0};
      clusters_M02 = {0};
      clusters_M20 = {0};
      clusters_isMatched = {0};
      clusters_NTower = {0};
      clusters_trueID = {0};
      clusters_NtrueID = {0};
      clusters_X = {0};
      clusters_Y = {0};
      clusters_Z = {0};
    }
  }
}


  const Int_t nBinsP              = 148;
  Double_t binningP[149]          = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  //10
                                      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,         // 20
                                      1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,         // 30
                                      2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,         // 40
                                      3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,         // 50
                                      4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0,         // 60
                                      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5,    // 70
                                      7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0,   // 80
                                      10.5, 11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15.0,   // 90
                                      15.5, 16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20.,    // 100
                                      21., 22., 23., 24., 25., 26., 27., 28., 29., 30.,         // 110
                                      31., 32., 33., 34., 35., 36., 37., 38., 39., 40.,         // 120
                                      42., 44., 46., 48., 50., 52., 54., 56., 58., 60.,         // 130
                                      65., 70., 75., 80., 85., 90., 95., 100., 110., 120.,      // 140
                                      130., 140., 150., 160., 170., 180., 190., 200., 250.};    // 149
                                       
