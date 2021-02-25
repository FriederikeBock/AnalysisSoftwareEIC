
const int _maxNHits = 5000;
const int _maxNTowers = 50*50;
const int _maxNclusters = 100;
const int _maxNProjections = 2000;
const int _maxNTracks = 200;
const int _maxNMCPart = 1000;


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
float* _tower_FHCAL_trueID       = new float[_maxNTowers];

// towers
int _nTowers_FEMC;
float* _tower_FEMC_E             = new float[_maxNTowers];
int* _tower_FEMC_iEta          = new int[_maxNTowers];
int* _tower_FEMC_iPhi          = new int[_maxNTowers];
float* _tower_FEMC_trueID        = new float[_maxNTowers];

// clusters
int _nclusters_FHCAL;
float* _cluster_FHCAL_E            = new float[_maxNclusters];
float* _cluster_FHCAL_Eta         = new float[_maxNclusters];
float* _cluster_FHCAL_Phi         = new float[_maxNclusters];
int* _cluster_FHCAL_NTower         = new int[_maxNclusters];
float* _cluster_FHCAL_trueID       = new float[_maxNclusters];

// clusters
int _nclusters_FEMC;
float* _cluster_FEMC_E             = new float[_maxNclusters];
float* _cluster_FEMC_Eta          = new float[_maxNclusters];
float* _cluster_FEMC_Phi          = new float[_maxNclusters];
int* _cluster_FEMC_NTower          = new int[_maxNclusters];
float* _cluster_FEMC_trueID        = new float[_maxNclusters];

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

    inputTree->SetBranchAddress("track_Proj_x",           _track_Proj_x);
    inputTree->SetBranchAddress("track_Proj_y",           _track_Proj_y);
    inputTree->SetBranchAddress("track_Proj_z",           _track_Proj_z);
    inputTree->SetBranchAddress("track_Proj_t",           _track_Proj_t);
    inputTree->SetBranchAddress("track_Proj_true_x",      _track_Proj_true_x);
    inputTree->SetBranchAddress("track_Proj_true_y",      _track_Proj_true_y);
    inputTree->SetBranchAddress("track_Proj_true_z",      _track_Proj_true_z);
    inputTree->SetBranchAddress("track_Proj_true_t",      _track_Proj_true_t);

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
    inputTree->SetBranchAddress("cluster_FHCAL_E",                _cluster_FHCAL_E);
    inputTree->SetBranchAddress("cluster_FHCAL_Eta",             _cluster_FHCAL_Eta);
    inputTree->SetBranchAddress("cluster_FHCAL_Phi",             _cluster_FHCAL_Phi);
    inputTree->SetBranchAddress("cluster_FHCAL_NTower",             _cluster_FHCAL_NTower);
    inputTree->SetBranchAddress("cluster_FHCAL_trueID",           _cluster_FHCAL_trueID);

    // clusters EMC
    inputTree->SetBranchAddress("cluster_FEMC_N",                 &_nclusters_FEMC);
    inputTree->SetBranchAddress("cluster_FEMC_E",                 _cluster_FEMC_E);
    inputTree->SetBranchAddress("cluster_FEMC_Eta",              _cluster_FEMC_Eta);
    inputTree->SetBranchAddress("cluster_FEMC_Phi",              _cluster_FEMC_Phi);
    inputTree->SetBranchAddress("cluster_FEMC_NTower",              _cluster_FEMC_NTower);
    inputTree->SetBranchAddress("cluster_FEMC_trueID",            _cluster_FEMC_trueID);

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