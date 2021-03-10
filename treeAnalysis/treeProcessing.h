TString outputDir;
const int _maxNHits = 5000;
const int _maxNTowers = 50*50;
const int _maxNclusters = 100;
const int _maxNProjections = 2000;
const int _maxNTracks = 200;
const int _maxNMCPart = 1000;

float _nEventsTree;


// track hits
int _nHitsLayers;
int* _hits_layerID              = new int[_maxNHits];
float* _hits_x             = new float[_maxNHits];
float* _hits_y             = new float[_maxNHits];
float* _hits_z             = new float[_maxNHits];
float* _hits_t             = new float[_maxNHits];

// towers
int _nTowers_FHCAL;
float* _tower_FHCAL_E            = new float[_maxNTowers];
int* _tower_FHCAL_iEta         = new int[_maxNTowers];
int* _tower_FHCAL_iPhi         = new int[_maxNTowers];
int* _tower_FHCAL_trueID       = new int[_maxNTowers];

// towers
int _nTowers_FEMC;
float* _tower_FEMC_E             = new float[_maxNTowers];
int* _tower_FEMC_iEta          = new int[_maxNTowers];
int* _tower_FEMC_iPhi          = new int[_maxNTowers];
int* _tower_FEMC_trueID        = new int[_maxNTowers];

// clusters
int _nclusters_FHCAL;
float* _clusters_FHCAL_E            = new float[_maxNclusters];
float* _clusters_FHCAL_Eta         = new float[_maxNclusters];
float* _clusters_FHCAL_Phi         = new float[_maxNclusters];
int* _clusters_FHCAL_NTower         = new int[_maxNclusters];
int* _clusters_FHCAL_trueID       = new int[_maxNclusters];

// clusters
int _nclusters_FEMC;
float* _clusters_FEMC_E             = new float[_maxNclusters];
float* _clusters_FEMC_Eta          = new float[_maxNclusters];
float* _clusters_FEMC_Phi          = new float[_maxNclusters];
int* _clusters_FEMC_NTower          = new int[_maxNclusters];
int* _clusters_FEMC_trueID        = new int[_maxNclusters];

// vertex
int _vertex_x;
int _vertex_y;
int _vertex_z;

// tracks
int _nTracks;
float* _track_ID                 = new float[_maxNTracks];
float* _track_trueID             = new float[_maxNTracks];
float* _track_px                 = new float[_maxNTracks];
float* _track_py                 = new float[_maxNTracks];
float* _track_pz                 = new float[_maxNTracks];


int _nProjections;
float* _track_ProjTrackID                 = new float[_maxNProjections];
int* _track_ProjLayer                 = new int[_maxNProjections];
float* _track_Proj_x             = new float[_maxNProjections];
float* _track_Proj_y             = new float[_maxNProjections];
float* _track_Proj_z             = new float[_maxNProjections];
float* _track_Proj_t             = new float[_maxNProjections];
float* _track_Proj_true_x             = new float[_maxNProjections];
float* _track_Proj_true_y             = new float[_maxNProjections];
float* _track_Proj_true_z             = new float[_maxNProjections];
float* _track_Proj_true_t             = new float[_maxNProjections];
// MC particles
int _nMCPart;
float* _mcpart_ID                = new float[_maxNMCPart];
float* _mcpart_ID_parent         = new float[_maxNMCPart];
float* _mcpart_PDG               = new float[_maxNMCPart];
float* _mcpart_E                 = new float[_maxNMCPart];
float* _mcpart_px                = new float[_maxNMCPart];
float* _mcpart_py                = new float[_maxNMCPart];
float* _mcpart_pz                = new float[_maxNMCPart];
float* _mcpart_Eta                = new float[_maxNMCPart];


void SetBranchAddressesTree(TTree* inputTree){
    inputTree->SetBranchAddress("nHits",                        &_nHitsLayers);
    inputTree->SetBranchAddress("hits_layerID",                 _hits_layerID);
    inputTree->SetBranchAddress("hits_x",               _hits_x);
    inputTree->SetBranchAddress("hits_y",               _hits_y);
    inputTree->SetBranchAddress("hits_z",               _hits_z);
    inputTree->SetBranchAddress("hits_t",               _hits_t);

    inputTree->SetBranchAddress("nTracks",                      &_nTracks);
    inputTree->SetBranchAddress("tracks_ID",                    _track_ID);
    inputTree->SetBranchAddress("tracks_px",                    _track_px);
    inputTree->SetBranchAddress("tracks_py",                    _track_py);
    inputTree->SetBranchAddress("tracks_pz",                    _track_pz);
    inputTree->SetBranchAddress("tracks_trueID",                _track_trueID);

    inputTree->SetBranchAddress("nProjections",        &_nProjections);
    inputTree->SetBranchAddress("track_ProjTrackID",   _track_ProjTrackID);
    inputTree->SetBranchAddress("track_ProjLayer",     _track_ProjLayer);

    inputTree->SetBranchAddress("track_TLP_x",           _track_Proj_x);
    inputTree->SetBranchAddress("track_TLP_y",           _track_Proj_y);
    inputTree->SetBranchAddress("track_TLP_z",           _track_Proj_z);
    inputTree->SetBranchAddress("track_TLP_t",           _track_Proj_t);
    inputTree->SetBranchAddress("track_TLP_true_x",      _track_Proj_true_x);
    inputTree->SetBranchAddress("track_TLP_true_y",      _track_Proj_true_y);
    inputTree->SetBranchAddress("track_TLP_true_z",      _track_Proj_true_z);
    inputTree->SetBranchAddress("track_TLP_true_t",      _track_Proj_true_t);

    // towers HCAL
    inputTree->SetBranchAddress("tower_FHCAL_N",                &_nTowers_FHCAL);
    inputTree->SetBranchAddress("tower_FHCAL_E",                _tower_FHCAL_E);
    inputTree->SetBranchAddress("tower_FHCAL_iEta",             _tower_FHCAL_iEta);
    inputTree->SetBranchAddress("tower_FHCAL_iPhi",             _tower_FHCAL_iPhi);
    inputTree->SetBranchAddress("tower_FHCAL_trueID",           _tower_FHCAL_trueID);

    // towers EMC
    inputTree->SetBranchAddress("tower_FEMC_N",                 &_nTowers_FEMC);
    inputTree->SetBranchAddress("tower_FEMC_E",                 _tower_FEMC_E);
    inputTree->SetBranchAddress("tower_FEMC_iEta",              _tower_FEMC_iEta);
    inputTree->SetBranchAddress("tower_FEMC_iPhi",              _tower_FEMC_iPhi);
    inputTree->SetBranchAddress("tower_FEMC_trueID",            _tower_FEMC_trueID);

    // clusters HCAL
    inputTree->SetBranchAddress("cluster_FHCAL_N",                &_nclusters_FHCAL);
    inputTree->SetBranchAddress("cluster_FHCAL_E",                _clusters_FHCAL_E);
    inputTree->SetBranchAddress("cluster_FHCAL_Eta",             _clusters_FHCAL_Eta);
    inputTree->SetBranchAddress("cluster_FHCAL_Phi",             _clusters_FHCAL_Phi);
    inputTree->SetBranchAddress("cluster_FHCAL_NTower",             _clusters_FHCAL_NTower);
    inputTree->SetBranchAddress("cluster_FHCAL_trueID",           _clusters_FHCAL_trueID);

    // clusters EMC
    inputTree->SetBranchAddress("cluster_FEMC_N",                 &_nclusters_FEMC);
    inputTree->SetBranchAddress("cluster_FEMC_E",                 _clusters_FEMC_E);
    inputTree->SetBranchAddress("cluster_FEMC_Eta",              _clusters_FEMC_Eta);
    inputTree->SetBranchAddress("cluster_FEMC_Phi",              _clusters_FEMC_Phi);
    inputTree->SetBranchAddress("cluster_FEMC_NTower",              _clusters_FEMC_NTower);
    inputTree->SetBranchAddress("cluster_FEMC_trueID",            _clusters_FEMC_trueID);

    // vertex
    inputTree->SetBranchAddress("vertex_x",                     &_vertex_x);
    inputTree->SetBranchAddress("vertex_y",                     &_vertex_y);
    inputTree->SetBranchAddress("vertex_z",                     &_vertex_z);

    // MC particles
    inputTree->SetBranchAddress("nMCPart",       &_nMCPart);
    inputTree->SetBranchAddress("mcpart_ID",     _mcpart_ID);
    inputTree->SetBranchAddress("mcpart_ID_parent",     _mcpart_ID_parent);
    inputTree->SetBranchAddress("mcpart_PDG",    _mcpart_PDG);
    inputTree->SetBranchAddress("mcpart_E",      _mcpart_E);
    inputTree->SetBranchAddress("mcpart_px",     _mcpart_px);
    inputTree->SetBranchAddress("mcpart_py",     _mcpart_py);
    inputTree->SetBranchAddress("mcpart_pz",     _mcpart_pz);
}

//__________________________________________________________________________________________________________
TString ReturnDateStr(){
    TDatime today;
    int iDate           = today.GetDate();
    int iYear           = iDate/10000;
    int iMonth          = (iDate%10000)/100;
    int iDay            = iDate%100;
    return Form("%i_%02d_%02d",iYear, iMonth, iDay);
}

const int _maxProjectionLayers = 5;
TString GetProjectionNameFromIndex(int projindex)
{
    switch (projindex)
    {
    case 0:    return "FTTL_0";
    case 1:    return "FTTL_1";
    case 2:    return "FTTL_2";
    case 3:    return "ETTL_0";
    case 4:    return "ETTL_1";
    // case 5:    return "FHCAL_0";
    // case 6:    return "FEMC_0";
    // case 7:    return "CTTL_0";
    // case 8:    return "CTTL_1";
    // case 9:    return "CTTL_2";
    default:   return "NOTHING";
    }
}