// ANCHOR basic struct for towers in clusterizer
typedef struct {
  float tower_E;
  int tower_iEta;
  int tower_iPhi;
  int tower_trueID;
} towersStrct;

// sorting function for towers
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.tower_E > rhs.tower_E; }


// conversion functions for processing
float* EtaPhiFromIndices(int ieta,int iphi, int granularity = 1,float energy = 0, bool isEMC = false);
TVector3 TowerPositionVectorFromIndices(int ieta,int iphi, bool isEMC = false, int granularity = 1);

float weightM02 = 4.5;
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, bool isEMC, Bool_t debugOutput);


// ANCHOR function to return eta and phi from cell iEta and iPhi indices
float * EtaPhiFromIndices(int i_x,int i_y, int granularity = 1, float energy, bool isEMC = false){
  static float eta_phi[2];
  float zHC = 400;
  float twrD = 10./granularity;
  // center (x=0,y=0) at 27,27 with 10x10 granularity
  int center_x = 27*granularity;
  int center_y = 27*granularity;
  if(isEMC){
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
TVector3 TowerPositionVectorFromIndices(int i_x,int i_y, bool isEMC = false, int granularity = 1){
  float zHC = 400;
  float twrD = 10./granularity;
  // center (x=0,y=0) at 27,27 with 10x10 granularity
  int center_x = 27*granularity;
  int center_y = 27*granularity;

  if(isEMC){
    zHC = 310;
    twrD = 5.535/granularity;
    center_x = 32*granularity;
    center_y = 32*granularity;
  }

  TVector3 twrPositionVec((center_x-i_x)*twrD,(center_y-i_y)*twrD,zHC);
  return twrPositionVec;
}


// ANCHOR function to determine shower shape
float * CalculateM02andWeightedPosition(std::vector<towersStrct> cluster_towers, float w_0, float cluster_E_calc, bool isEMC, Bool_t debugOutput){
    static float returnVariables[4]; //0:M02, 1:M20, 2:eta, 3: phi
    float w_tot = 0;
    std::vector<float> w_i;
    TVector3 vecTwr;

    //calculation of weights and weighted position vector
    for(int cellI=0; cellI<cluster_towers.size(); cellI++){
        w_i.push_back(TMath::Max( (float)0, (float) (w_0 + TMath::Log(cluster_towers.at(cellI).tower_E/cluster_E_calc) )));
        w_tot += w_i.at(cellI);
        vecTwr += w_i.at(cellI)*TowerPositionVectorFromIndices(cluster_towers.at(cellI).tower_iEta,cluster_towers.at(cellI).tower_iPhi, isEMC);
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

