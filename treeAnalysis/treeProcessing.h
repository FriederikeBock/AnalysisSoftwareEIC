#include <TROOT.h>
#include <TString.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

TString outputDir;
// Constants are from the EventEvaluator class in coresoftware
const int _maxNHits = 10000;
const int _maxNTowers = 50 * 50;
const int _maxNTowersCentral = 2000;
const int _maxNTowersDR = 3000 * 3000;
const int _maxNTowersCalo = 5000000;
const int _maxNclusters = 100;
const int _maxNclustersCentral = 2000;
const int _maxNTracks = 200;
const int _maxNProjections = 2000;
const int _maxNMCPart = 100000;
const int _maxNHepmcp = 1000;
int verbosityBASE = 2;

float _nEventsTree;

enum calotype {
  kFHCAL          = 0,
  kFEMC           = 1,
  kDRCALO         = 2,
  kEEMC           = 3,
  kCEMC           = 4,
  kEHCAL          = 5,
  kHCALIN         = 6,
  kHCALOUT        = 7,
  kLFHCAL         = 8,
  kEEMCG          = 9,
  kBECAL          = 10,
  kFOCAL          = 11
};

bool caloEnabled[20]      = {0};
bool tracksEnabled        = 0;
bool vertexEnabled        = 0;
bool xSectionEnabled      = 0;

// Event level info
float _cross_section;
float _event_weight;
int _n_generator_accepted;

// track hits
int _nHitsLayers;
int* _hits_layerID              = new int[_maxNHits];
float* _hits_x             = new float[_maxNHits];
float* _hits_y             = new float[_maxNHits];
float* _hits_z             = new float[_maxNHits];
float* _hits_t             = new float[_maxNHits];
float* _hits_edep          = new float[_maxNHits];

// towers
int _nTowers_CEMC;
float* _tower_CEMC_E            = new float[_maxNTowersCentral];
int* _tower_CEMC_iEta         = new int[_maxNTowersCentral];
int* _tower_CEMC_iPhi         = new int[_maxNTowersCentral];
int* _tower_CEMC_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EEMC;
float* _tower_EEMC_E            = new float[_maxNTowersCentral];
int* _tower_EEMC_iEta         = new int[_maxNTowersCentral];
int* _tower_EEMC_iPhi         = new int[_maxNTowersCentral];
int* _tower_EEMC_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EEMCG;
float* _tower_EEMCG_E            = new float[_maxNTowersCentral];
int* _tower_EEMCG_iEta         = new int[_maxNTowersCentral];
int* _tower_EEMCG_iPhi         = new int[_maxNTowersCentral];
int* _tower_EEMCG_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EHCAL;
float* _tower_EHCAL_E            = new float[_maxNTowersCentral];
int* _tower_EHCAL_iEta         = new int[_maxNTowersCentral];
int* _tower_EHCAL_iPhi         = new int[_maxNTowersCentral];
int* _tower_EHCAL_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALIN;
float* _tower_HCALIN_E            = new float[_maxNTowersCentral];
int* _tower_HCALIN_iEta         = new int[_maxNTowersCentral];
int* _tower_HCALIN_iPhi         = new int[_maxNTowersCentral];
int* _tower_HCALIN_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALOUT;
float* _tower_HCALOUT_E            = new float[_maxNTowersCentral];
int* _tower_HCALOUT_iEta         = new int[_maxNTowersCentral];
int* _tower_HCALOUT_iPhi         = new int[_maxNTowersCentral];
int* _tower_HCALOUT_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_DRCALO;
float* _tower_DRCALO_E            = new float[_maxNTowersDR];
int* _tower_DRCALO_iEta         = new int[_maxNTowersDR];
int* _tower_DRCALO_iPhi         = new int[_maxNTowersDR];
int* _tower_DRCALO_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_FOCAL;
float* _tower_FOCAL_E            = new float[_maxNTowersDR];
int* _tower_FOCAL_iEta         = new int[_maxNTowersDR];
int* _tower_FOCAL_iPhi         = new int[_maxNTowersDR];
int* _tower_FOCAL_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_LFHCAL;
float* _tower_LFHCAL_E            = new float[_maxNTowersDR];
int* _tower_LFHCAL_iEta         = new int[_maxNTowersDR];
int* _tower_LFHCAL_iPhi         = new int[_maxNTowersDR];
int* _tower_LFHCAL_iL           = new int[_maxNTowersDR];
int* _tower_LFHCAL_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_BECAL;
float* _tower_BECAL_E            = new float[_maxNTowersDR];
int* _tower_BECAL_iEta         = new int[_maxNTowersDR];
int* _tower_BECAL_iPhi         = new int[_maxNTowersDR];
int* _tower_BECAL_trueID       = new int[_maxNTowersDR];

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
float _vertex_x;
float _vertex_y;
float _vertex_z;

// tracks
int _nTracks;
float* _track_ID                 = new float[_maxNTracks];
float* _track_trueID             = new float[_maxNTracks];
float* _track_px                 = new float[_maxNTracks];
float* _track_py                 = new float[_maxNTracks];
float* _track_pz                 = new float[_maxNTracks];
unsigned short* _track_source             = new unsigned short[_maxNTracks];
std::array<std::vector<int>, _maxNTracks> _track_RefProjID;
// Initializes all elements to 0
std::array<bool, _maxNTracks> _track_hasTTL{{}};
std::array<int, _maxNTracks> _track_nTTL{{}};
std::array<int, _maxNTracks> _track_nTrL{{}};
std::array<bool, _maxNTracks> _track_hasIL{{}};
std::array<bool, _maxNTracks> _track_hasOL{{}};

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
std::array<int, _maxNProjections> _track_Proj_Clas{{}};
// MC particles
int _nMCPart;
int* _mcpart_ID                   = new int[_maxNMCPart];
int* _mcpart_ID_parent            = new int[_maxNMCPart];
int* _mcpart_PDG                  = new int[_maxNMCPart];
float* _mcpart_E                  = new float[_maxNMCPart];
float* _mcpart_px                 = new float[_maxNMCPart];
float* _mcpart_py                 = new float[_maxNMCPart];
float* _mcpart_pz                 = new float[_maxNMCPart];
float* _mcpart_Eta                = new float[_maxNMCPart];
float* _mcpart_Phi                = new float[_maxNMCPart];
std::array<std::vector<int>, _maxNMCPart> _mcpart_RecTrackIDs;


TRandom3  _fRandom;                                  // random for effi generation


int _calogeom_ID;
int _calogeom_towers_N;
int*  _calogeom_towers_iEta = new int[_maxNTowersCalo];
int*  _calogeom_towers_iPhi = new int[_maxNTowersCalo];
int*  _calogeom_towers_iL   = new int[_maxNTowersCalo];
float*  _calogeom_towers_Eta = new float[_maxNTowersCalo];
float*  _calogeom_towers_Phi = new float[_maxNTowersCalo];
float*  _calogeom_towers_x = new float[_maxNTowersCalo];
float*  _calogeom_towers_y = new float[_maxNTowersCalo];
float*  _calogeom_towers_z = new float[_maxNTowersCalo];

TTree * tt_geometry;

void SetBranchAddressesGeometryTree(TTree* inputTreeGeo){
    inputTreeGeo->SetBranchAddress("calo",              &_calogeom_ID);
    inputTreeGeo->SetBranchAddress("calo_towers_N",     &_calogeom_towers_N);
    inputTreeGeo->SetBranchAddress("calo_towers_iEta",  _calogeom_towers_iEta);
    inputTreeGeo->SetBranchAddress("calo_towers_iPhi",  _calogeom_towers_iPhi);
    inputTreeGeo->SetBranchAddress("calo_towers_iL",    _calogeom_towers_iL);
    inputTreeGeo->SetBranchAddress("calo_towers_Eta",   _calogeom_towers_Eta);
    inputTreeGeo->SetBranchAddress("calo_towers_Phi",   _calogeom_towers_Phi);
    inputTreeGeo->SetBranchAddress("calo_towers_x",     _calogeom_towers_x);
    inputTreeGeo->SetBranchAddress("calo_towers_y",     _calogeom_towers_y);
    inputTreeGeo->SetBranchAddress("calo_towers_z",     _calogeom_towers_z);
}
void SetBranchAddressesTree(TTree* inputTree){

    if (inputTree->GetBranchStatus("cross_section") ){
      xSectionEnabled = 1;
      inputTree->SetBranchAddress("cross_section", &_cross_section);
      inputTree->SetBranchAddress("event_weight", &_event_weight);
      inputTree->SetBranchAddress("n_generator_accepted", &_n_generator_accepted);
    }
    if (inputTree->GetBranchStatus("nHits") ){
      inputTree->SetBranchAddress("nHits",                        &_nHitsLayers);
      inputTree->SetBranchAddress("hits_layerID",                 _hits_layerID);
      inputTree->SetBranchAddress("hits_x",               _hits_x);
      inputTree->SetBranchAddress("hits_y",               _hits_y);
      inputTree->SetBranchAddress("hits_z",               _hits_z);
      inputTree->SetBranchAddress("hits_t",               _hits_t);
      inputTree->SetBranchAddress("hits_edep",               _hits_edep);
    }
    
    if (inputTree->GetBranchStatus("nTracks") ){
      tracksEnabled = 1;
      inputTree->SetBranchAddress("nTracks",              &_nTracks);
      inputTree->SetBranchAddress("tracks_ID",            _track_ID);
      inputTree->SetBranchAddress("tracks_px",            _track_px);
      inputTree->SetBranchAddress("tracks_py",            _track_py);
      inputTree->SetBranchAddress("tracks_pz",            _track_pz);
      inputTree->SetBranchAddress("tracks_trueID",        _track_trueID);
      inputTree->SetBranchAddress("tracks_source",        _track_source);
      
      inputTree->SetBranchAddress("nProjections",         &_nProjections);
      inputTree->SetBranchAddress("track_ProjTrackID",    _track_ProjTrackID);
      inputTree->SetBranchAddress("track_ProjLayer",      _track_ProjLayer);

      inputTree->SetBranchAddress("track_TLP_x",           _track_Proj_x);
      inputTree->SetBranchAddress("track_TLP_y",           _track_Proj_y);
      inputTree->SetBranchAddress("track_TLP_z",           _track_Proj_z);
      inputTree->SetBranchAddress("track_TLP_t",           _track_Proj_t);
      inputTree->SetBranchAddress("track_TLP_true_x",      _track_Proj_true_x);
      inputTree->SetBranchAddress("track_TLP_true_y",      _track_Proj_true_y);
      inputTree->SetBranchAddress("track_TLP_true_z",      _track_Proj_true_z);
      inputTree->SetBranchAddress("track_TLP_true_t",      _track_Proj_true_t);
    }

    // towers EEMC
    if( inputTree->GetBranchStatus("tower_EEMC_N") ){
      caloEnabled[kEEMC] = 1;
      inputTree->SetBranchAddress("tower_EEMC_N",                &_nTowers_EEMC);
      inputTree->SetBranchAddress("tower_EEMC_E",                _tower_EEMC_E);
      inputTree->SetBranchAddress("tower_EEMC_iEta",             _tower_EEMC_iEta);
      inputTree->SetBranchAddress("tower_EEMC_iPhi",             _tower_EEMC_iPhi);
      inputTree->SetBranchAddress("tower_EEMC_trueID",           _tower_EEMC_trueID);
    } 
    // towers EEMCG
    if( inputTree->GetBranchStatus("tower_EEMCG_N") ){
      caloEnabled[kEEMCG] = 1;
      inputTree->SetBranchAddress("tower_EEMCG_N",                &_nTowers_EEMCG);
      inputTree->SetBranchAddress("tower_EEMCG_E",                _tower_EEMCG_E);
      inputTree->SetBranchAddress("tower_EEMCG_iEta",             _tower_EEMCG_iEta);
      inputTree->SetBranchAddress("tower_EEMCG_iPhi",             _tower_EEMCG_iPhi);
      inputTree->SetBranchAddress("tower_EEMCG_trueID",           _tower_EEMCG_trueID);
    } 

    // towers EHCAL
    if( inputTree->GetBranchStatus("tower_EHCAL_N") ){
      caloEnabled[kEHCAL] = 1;
      inputTree->SetBranchAddress("tower_EHCAL_N",                &_nTowers_EHCAL);
      inputTree->SetBranchAddress("tower_EHCAL_E",                _tower_EHCAL_E);
      inputTree->SetBranchAddress("tower_EHCAL_iEta",             _tower_EHCAL_iEta);
      inputTree->SetBranchAddress("tower_EHCAL_iPhi",             _tower_EHCAL_iPhi);
      inputTree->SetBranchAddress("tower_EHCAL_trueID",           _tower_EHCAL_trueID);
    }

    // towers HCALIN
    if( inputTree->GetBranchStatus("tower_HCALIN_N") ){
      caloEnabled[kHCALIN] = 1;
      inputTree->SetBranchAddress("tower_HCALIN_N",                &_nTowers_HCALIN);
      inputTree->SetBranchAddress("tower_HCALIN_E",                _tower_HCALIN_E);
      inputTree->SetBranchAddress("tower_HCALIN_iEta",             _tower_HCALIN_iEta);
      inputTree->SetBranchAddress("tower_HCALIN_iPhi",             _tower_HCALIN_iPhi);
      inputTree->SetBranchAddress("tower_HCALIN_trueID",           _tower_HCALIN_trueID);
    } 

    // towers HCALOUT
    if( inputTree->GetBranchStatus("tower_HCALOUT_N") ){
      caloEnabled[kHCALOUT] = 1;
      inputTree->SetBranchAddress("tower_HCALOUT_N",                &_nTowers_HCALOUT);
      inputTree->SetBranchAddress("tower_HCALOUT_E",                _tower_HCALOUT_E);
      inputTree->SetBranchAddress("tower_HCALOUT_iEta",             _tower_HCALOUT_iEta);
      inputTree->SetBranchAddress("tower_HCALOUT_iPhi",             _tower_HCALOUT_iPhi);
      inputTree->SetBranchAddress("tower_HCALOUT_trueID",           _tower_HCALOUT_trueID);
    } 

    // towers CEMC
    if( inputTree->GetBranchStatus("tower_CEMC_N") ){
      caloEnabled[kCEMC] = 1;
      inputTree->SetBranchAddress("tower_CEMC_N",                &_nTowers_CEMC);
      inputTree->SetBranchAddress("tower_CEMC_E",                _tower_CEMC_E);
      inputTree->SetBranchAddress("tower_CEMC_iEta",             _tower_CEMC_iEta);
      inputTree->SetBranchAddress("tower_CEMC_iPhi",             _tower_CEMC_iPhi);
      inputTree->SetBranchAddress("tower_CEMC_trueID",           _tower_CEMC_trueID);
    } 

    // towers BECAL
    if( inputTree->GetBranchStatus("tower_BECAL_N") ){
      caloEnabled[kBECAL] = 1;
      inputTree->SetBranchAddress("tower_BECAL_N",                &_nTowers_BECAL);
      inputTree->SetBranchAddress("tower_BECAL_E",                _tower_BECAL_E);
      inputTree->SetBranchAddress("tower_BECAL_iEta",             _tower_BECAL_iEta);
      inputTree->SetBranchAddress("tower_BECAL_iPhi",             _tower_BECAL_iPhi);
      inputTree->SetBranchAddress("tower_BECAL_trueID",           _tower_BECAL_trueID);
    } 

    // towers LFHCAL
    if( inputTree->GetBranchStatus("tower_LFHCAL_N") ){
      caloEnabled[kLFHCAL] = 1;
      inputTree->SetBranchAddress("tower_LFHCAL_N",                &_nTowers_LFHCAL);
      inputTree->SetBranchAddress("tower_LFHCAL_E",                _tower_LFHCAL_E);
      inputTree->SetBranchAddress("tower_LFHCAL_iEta",             _tower_LFHCAL_iEta);
      inputTree->SetBranchAddress("tower_LFHCAL_iPhi",             _tower_LFHCAL_iPhi);
      inputTree->SetBranchAddress("tower_LFHCAL_iL",               _tower_LFHCAL_iL);
      inputTree->SetBranchAddress("tower_LFHCAL_trueID",           _tower_LFHCAL_trueID);
    } 
    // towers DRCALO
    if( inputTree->GetBranchStatus("tower_DRCALO_N") ){
      caloEnabled[kDRCALO] = 1;
      inputTree->SetBranchAddress("tower_DRCALO_N",                &_nTowers_DRCALO);
      inputTree->SetBranchAddress("tower_DRCALO_E",                _tower_DRCALO_E);
      inputTree->SetBranchAddress("tower_DRCALO_iEta",             _tower_DRCALO_iEta);
      inputTree->SetBranchAddress("tower_DRCALO_iPhi",             _tower_DRCALO_iPhi);
      inputTree->SetBranchAddress("tower_DRCALO_trueID",           _tower_DRCALO_trueID);
    } 
    // towers DRCALO
    if( inputTree->GetBranchStatus("tower_FOCAL_N") ){
      caloEnabled[kFOCAL] = 1;
      inputTree->SetBranchAddress("tower_FOCAL_N",                &_nTowers_FOCAL);
      inputTree->SetBranchAddress("tower_FOCAL_E",                _tower_FOCAL_E);
      inputTree->SetBranchAddress("tower_FOCAL_iEta",             _tower_FOCAL_iEta);
      inputTree->SetBranchAddress("tower_FOCAL_iPhi",             _tower_FOCAL_iPhi);
      inputTree->SetBranchAddress("tower_FOCAL_trueID",           _tower_FOCAL_trueID);
    } 

    // towers FHCAL
    if( inputTree->GetBranchStatus("tower_FHCAL_N") ){
      caloEnabled[kFHCAL] = 1;
      inputTree->SetBranchAddress("tower_FHCAL_N",                &_nTowers_FHCAL);
      inputTree->SetBranchAddress("tower_FHCAL_E",                _tower_FHCAL_E);
      inputTree->SetBranchAddress("tower_FHCAL_iEta",             _tower_FHCAL_iEta);
      inputTree->SetBranchAddress("tower_FHCAL_iPhi",             _tower_FHCAL_iPhi);
      inputTree->SetBranchAddress("tower_FHCAL_trueID",           _tower_FHCAL_trueID);
    } 

    // towers FEMC
    if( inputTree->GetBranchStatus("tower_FEMC_N") ){
      caloEnabled[kFEMC] = 1;
      inputTree->SetBranchAddress("tower_FEMC_N",                 &_nTowers_FEMC);
      inputTree->SetBranchAddress("tower_FEMC_E",                 _tower_FEMC_E);
      inputTree->SetBranchAddress("tower_FEMC_iEta",              _tower_FEMC_iEta);
      inputTree->SetBranchAddress("tower_FEMC_iPhi",              _tower_FEMC_iPhi);
      inputTree->SetBranchAddress("tower_FEMC_trueID",            _tower_FEMC_trueID);
    }

    // clusters HCAL
    if (caloEnabled[kFHCAL]){
      if (inputTree->GetBranchStatus("cluster_FHCAL_N") ){
        inputTree->SetBranchAddress("cluster_FHCAL_N",                &_nclusters_FHCAL);
        inputTree->SetBranchAddress("cluster_FHCAL_E",                _clusters_FHCAL_E);
        inputTree->SetBranchAddress("cluster_FHCAL_Eta",             _clusters_FHCAL_Eta);
        inputTree->SetBranchAddress("cluster_FHCAL_Phi",             _clusters_FHCAL_Phi);
        inputTree->SetBranchAddress("cluster_FHCAL_NTower",             _clusters_FHCAL_NTower);
        inputTree->SetBranchAddress("cluster_FHCAL_trueID",           _clusters_FHCAL_trueID);
      }
    }
    // clusters EMC
    if (caloEnabled[kFEMC]){
      if (inputTree->GetBranchStatus("cluster_FEMC_N") ){
        inputTree->SetBranchAddress("cluster_FEMC_N",                 &_nclusters_FEMC);
        inputTree->SetBranchAddress("cluster_FEMC_E",                 _clusters_FEMC_E);
        inputTree->SetBranchAddress("cluster_FEMC_Eta",              _clusters_FEMC_Eta);
        inputTree->SetBranchAddress("cluster_FEMC_Phi",              _clusters_FEMC_Phi);
        inputTree->SetBranchAddress("cluster_FEMC_NTower",              _clusters_FEMC_NTower);
        inputTree->SetBranchAddress("cluster_FEMC_trueID",            _clusters_FEMC_trueID);
      }
    }
    // vertex

    if (inputTree->GetBranchStatus("vertex_x") ){
      vertexEnabled = 1;
      inputTree->SetBranchAddress("vertex_x",                     &_vertex_x);
      inputTree->SetBranchAddress("vertex_y",                     &_vertex_y);
      inputTree->SetBranchAddress("vertex_z",                     &_vertex_z);
    }
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

bool _do_TimingStudies        = false;
bool _is_ALLSILICON           = false;
const int _maxProjectionLayers    = 32;
Int_t layerIndexHist[32]  =  { 0, 1, 2, 3, 4, 7,        // TTLs
                                 50, 51, 52, 53, 54,      // FST
                                 100, 101, 102, 103,      // EFST
                                 155, 156, 157, 150, 151, // SVTX & BARR
                                 110, 111, 112, 120, 130,  // mu RWell 
                                 6, 60, 61, 62, 63, 66, 67 // calo's 
};

// void ResetLayerIndexForward(){
//   if (!_is_ALLSILICON){
//     layerIndexHist[0]  = 50;
//     layerIndexHist[1]  = 51;
//     layerIndexHist[2]  = 52;
//     layerIndexHist[3]  = 53;
//     layerIndexHist[4]  = 54;
//   }
// }

TString GetProjectionNameFromIndex(int projindex)
{
  switch (projindex){
    // timing layers
    case 0:     return "FTTL_0";
    case 1:     return "FTTL_1";
    case 2:     return "FTTL_2";
    case 3:     return "ETTL_0";
    case 4:     return "ETTL_1";
    case 7:     return "CTTL_0";
    case 8:     return "CTTL_1";
    // LBL central barrel
    case 10:    return "CENTAL_0"; //"LBLVTX_CENTRAL_10";
    case 11:    return "CENTAL_1"; //"LBLVTX_CENTRAL_11";
    case 12:    return "CENTAL_2"; //"LBLVTX_CENTRAL_12";
    case 13:    return "CENTAL_3"; //"LBLVTX_CENTRAL_13";
    case 14:    return "CENTAL_4"; //"LBLVTX_CENTRAL_14";
    case 15:    return "CENTAL_5"; //"LBLVTX_CENTRAL_15";
    // LBL forward barrel
    case 20:    return "FORWARD_0"; // "LBLVTX_FORWARD_20"
    case 21:    return "FORWARD_1"; // "LBLVTX_FORWARD_21"
    case 22:    return "FORWARD_2"; // "LBLVTX_FORWARD_22"
    case 23:    return "FORWARD_3"; // "LBLVTX_FORWARD_23"
    case 24:    return "FORWARD_4"; // "LBLVTX_FORWARD_24"
    // LBL backward barrel
    case 30:    return "BACKWARD_0"; //"LBLVTX_BACKWARD_30";
    case 31:    return "BACKWARD_1"; //"LBLVTX_BACKWARD_31";
    case 32:    return "BACKWARD_2"; //"LBLVTX_BACKWARD_32";
    case 33:    return "BACKWARD_3"; //"LBLVTX_BACKWARD_33";
    case 34:    return "BACKWARD_4"; //"LBLVTX_BACKWARD_34";
  
    // old BARREL tracker
    case 40:    return "BARREL_0"; //"BARREL_0";
    case 41:    return "BARREL_1"; //"BARREL_1";
    case 42:    return "BARREL_2"; //"BARREL_2";
    case 43:    return "BARREL_3"; //"BARREL_3";
    case 44:    return "BARREL_4"; //"BARREL_4";
    case 45:    return "BARREL_5"; //"BARREL_5";
    // FST 
    case 50:    return "FST_0";
    case 51:    return "FST_1";
    case 52:    return "FST_2";
    case 53:    return "FST_3";
    case 54:    return "FST_4";
    case 55:    return "FST_5";
    // EFST
    case 100:    return "EFST_0";
    case 101:    return "EFST_1";
    case 102:    return "EFST_2";
    case 103:    return "EFST_3";
    case 104:    return "EFST_4";
    case 105:    return "EFST_5";
    // new BARREL tracker
    case 150:    return "BARR_0"; 
    case 151:    return "BARR_1"; 
    case 152:    return "BARR_2"; 
    case 153:    return "BARR_3"; 
    case 154:    return "BARR_4"; 
    case 160:    return "BARR"; 
    // SVTX
    case 155:    return "SVTX_0"; 
    case 156:    return "SVTX_1"; 
    case 157:    return "SVTX_2"; 
    case 158:    return "SVTX_3"; 
    case 159:    return "SVTX_4"; 
    case 161:    return "SVTX"; 
    // mu RWells
    case 110:    return "RWELL_0"; 
    case 111:    return "RWELL_1"; 
    case 112:    return "RWELL_2"; 
    
    // forward GEMS
    case 120:    return "FGEM_0"; 
    case 121:    return "FGEM_1"; 
    // backward GEMS
    case 130:    return "EGEM_0"; 
    case 131:    return "EGEM_1"; 
    
    // calorimeters
    case 5:    return "FHCAL";
    case 6:    return "FEMC";
    case 60:   return "EHCAL";
    case 61:   return "EEMC";
    case 62:   return "HCALIN";
    case 63:   return "HCALOUT";
    case 64:   return "CEMC";
    case 65:   return "EEMC_glass_0";
    case 66:   return "BECAL";
    case 67:   return "LFHCAL";
    case 140:  return "LFHCAL_0";
    case 141:  return "LFHCAL_1";
    case 142:  return "LFHCAL_2";
    case 143:  return "LFHCAL_3";
    case 144:  return "LFHCAL_4";
    case 145:  return "LFHCAL_5";
    case 146:  return "LFHCAL_6";
    case 147:  return "LFHCAL_7";
    
    default:   return "NOTHING";
  }
}

Int_t GetRegionFromIndex(int projindex)
{
  switch (projindex){
    // timing layers
    case 0:     return 2; // "FTTL_0";
    case 1:     return 2; //"FTTL_1";
    case 2:     return 2; //"FTTL_2";
    case 3:     return 0; // "ETTL_0";
    case 4:     return 0; // "ETTL_1";
    case 7:     return 1; // "CTTL_0";
    case 8:     return 1; // "CTTL_1";
    // LBL central barrel
    case 10:    return 1; // "CENTAL_0"; //"LBLVTX_CENTRAL_10";
    case 11:    return 1; // "CENTAL_1"; //"LBLVTX_CENTRAL_11";
    case 12:    return 1; // "CENTAL_2"; //"LBLVTX_CENTRAL_12";
    case 13:    return 1; // "CENTAL_3"; //"LBLVTX_CENTRAL_13";
    case 14:    return 1; // "CENTAL_4"; //"LBLVTX_CENTRAL_14";
    case 15:    return 1; // "CENTAL_5"; //"LBLVTX_CENTRAL_15";
    // LBL forward barrel
    case 20:    return 2; //"FORWARD_0"; // "LBLVTX_FORWARD_20"
    case 21:    return 2; //"FORWARD_1"; // "LBLVTX_FORWARD_21"
    case 22:    return 2; //"FORWARD_2"; // "LBLVTX_FORWARD_22"
    case 23:    return 2; //"FORWARD_3"; // "LBLVTX_FORWARD_23"
    case 24:    return 2; //"FORWARD_4"; // "LBLVTX_FORWARD_24"
    // LBL backward barrel
    case 30:    return 0; //"BACKWARD_0"; //"LBLVTX_BACKWARD_30";
    case 31:    return 0; //"BACKWARD_1"; //"LBLVTX_BACKWARD_31";
    case 32:    return 0; //"BACKWARD_2"; //"LBLVTX_BACKWARD_32";
    case 33:    return 0; //"BACKWARD_3"; //"LBLVTX_BACKWARD_33";
    case 34:    return 0; //"BACKWARD_4"; //"LBLVTX_BACKWARD_34";
  
    // old BARREL tracker
    case 40:    return 1; // "BARREL_0"; //"BARREL_0";
    case 41:    return 1; // "BARREL_1"; //"BARREL_1";
    case 42:    return 1; // "BARREL_2"; //"BARREL_2";
    case 43:    return 1; // "BARREL_3"; //"BARREL_3";
    case 44:    return 1; // "BARREL_4"; //"BARREL_4";
    case 45:    return 1; // "BARREL_5"; //"BARREL_5";
    // FST 
    case 50:    return 2; // "FST_0";
    case 51:    return 2; // "FST_1";
    case 52:    return 2; // "FST_2";
    case 53:    return 2; // "FST_3";
    case 54:    return 2; // "FST_4";
    case 55:    return 2; // "FST_5";
    // EFST
    case 100:    return 0; // "EFST_0";
    case 101:    return 0; // "EFST_1";
    case 102:    return 0; // "EFST_2";
    case 103:    return 0; // "EFST_3";
    case 104:    return 0; // "EFST_4";
    case 105:    return 0; // "EFST_5";
    // new BARREL tracker
    case 150:    return 1; // "BARR_0"; 
    case 151:    return 1; // "BARR_1"; 
    case 152:    return 1; // "BARR_2"; 
    case 153:    return 1; // "BARR_3"; 
    case 154:    return 1; // "BARR_4"; 
    case 160:    return 1; // "BARR"; 
    // SVTX
    case 155:    return 1; // "SVTX_0"; 
    case 156:    return 1; // "SVTX_1"; 
    case 157:    return 1; // "SVTX_2"; 
    case 158:    return 1; // "SVTX_3"; 
    case 159:    return 1; // "SVTX_4"; 
    case 161:    return 1; // "SVTX"; 
    // mu RWells
    case 110:    return 1; // "RWELL_1"; 
    case 111:    return 1; // "RWELL_1"; 
    case 112:    return 1; // "RWELL_1"; 
    
    // forward GEMS
    case 120:    return 2; // "FGEM_0"; 
    case 121:    return 2; // "FGEM_1"; 
    // backward GEMS
    case 130:    return 0; // "EGEM_0"; 
    case 131:    return 0; // "EGEM_1"; 
    
    // calorimeters
    case 5:    return 2; // "FHCAL";
    case 6:    return 2; // "FEMC";
    case 60:   return 0; // "EHCAL";
    case 61:   return 0; // "EEMC";
    case 62:   return 1; // "HCALIN";
    case 63:   return 1; // "HCALOUT";
    case 64:   return 1; // "CEMC";
    case 65:   return 0; // "EEMC_glass_0";
    case 66:   return 1; // "BECAL";
    case 67:   return 2; // "LFHCAL";
    case 140:  return 2; // "LFHCAL_0";
    case 141:  return 2; // "LFHCAL_1";
    case 142:  return 2; // "LFHCAL_2";
    case 143:  return 2; // "LFHCAL_3";
    case 144:  return 2; // "LFHCAL_4";
    case 145:  return 2; // "LFHCAL_5";
    case 146:  return 2; // "LFHCAL_6";
    case 147:  return 2; // "LFHCAL_7";
    
    default:   return 0;
  }
}

Int_t ReturnIndexForwardLayer(Int_t layerID){
  for (Int_t i = 0; i < _maxProjectionLayers; i++){
    if (layerIndexHist[i] == layerID)
      return i;
  }
  return -1;
}


Bool_t HasFirstTwoLayers(Int_t layerID){
  if (_is_ALLSILICON){
    switch (layerID){  
      case 20: 
      case 21: 
      case 30:
      case 31:
      case 10:
      case 11:
        return kTRUE;
      default:
        return kFALSE;
    }
  } else {
    switch (layerID){  
      case 50:  // FST_0
      case 51:  // FST_1
      case 100: // EFST_0
      case 101: // EFST_1
      case 154: // SVTX_0
      case 155: // SVTX_1
      case 156: // SVTX_2
      case 161: // SVTX
      case 40:  // BARREL_0
      case 41:  // BARREL_1
        return kTRUE;
      default:
        return kFALSE;
    }    
  }
  return kFALSE;
}

Bool_t IsTrackerLayer(Int_t layerID){
  if (_is_ALLSILICON){
    switch (layerID){  
      case 20: 
      case 21: 
      case 22: 
      case 23: 
      case 24: 
      case 30:
      case 31:
      case 32:
      case 33:
      case 34:
      case 10:
      case 11:
      case 12:
      case 13:
      case 14:
      case 15:
        return kTRUE;
      default:
        return kFALSE;
    }
  } else {
    switch (layerID){  
      // FST
      case 50: 
      case 51: 
      case 52: 
      case 53: 
      case 54: 
      case 55:   
      // BARR by layer
      case 150:
      case 151:
      case 152:
      case 153:
      case 154:
      // SVTX by layer
      case 155:
      case 156:
      case 157:
      case 158:
      case 159:
      // BARR
      case 160:
      // SVTX
      case 161:
      // EFST
      case 100:
      case 101:
      case 102:
      case 103:
      case 104:
      case 105:
      case 106:
      // RWELL
      case 110:
      case 111:
      case 112:
      // FGEM
      case 120:
      case 121:
      // EGEM
      case 130:
      case 131:
        return kTRUE;
      default:
        return kFALSE;
    }    
  }
  return kFALSE;
}

Bool_t IsCaloProjection(Int_t layerID){
  switch (layerID){  
    case 5: 
    case 6: 
    case 60: 
    case 61: 
    case 62: 
    case 63:
    case 64:
    case 65:
    case 66:
    case 67:
    case 140:
    case 141:
    case 142:
    case 143:
    case 144:
    case 145:
    case 146:
    case 147:
      return kTRUE;
    default:
      return kFALSE;
  }
  return kFALSE;
}

Bool_t IsFarForwardProjection(Int_t layerID){
  switch (layerID){  
    case 70: 
    case 71: 
    case 72: 
    case 73: 
    case 74: 
    case 90:
    case 91:
    case 92:
      return kTRUE;
    default:
      return kFALSE;
  }
  return kFALSE;
}

Bool_t HasTimingLayer(Int_t layerID){
  switch (layerID){
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 7:
    case 8:
      return kTRUE;
    default:
      return kFALSE;
  }
  return kFALSE;
}

// **********************************************************************************************
// ****************** resetting of MC references  ***********************************************
// **********************************************************************************************
int GetCorrectMCArrayEntry(float objectTrueID){
  for(Int_t imc=0; imc<_nMCPart; imc++){
    if(objectTrueID==_mcpart_ID[imc]){
      return imc;
    }
  }
  return -1;
}


// **********************************************************************************************
// ****************** create vectors for matching diff rec tracks to MC part ********************
// **********************************************************************************************
void prepareMCMatchInfo(){
  for(Int_t itrk=0; itrk<(Int_t)_nTracks; itrk++){
    if (verbosityBASE > 3) std::cout << "processing track: " << itrk << std::endl;
    if (_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size() > 0){
      if (verbosityBASE > 3)
        std::cout << "adding rec track: \t" << itrk << "\t track source: " << _track_source[itrk] << " MC: "<< (int)_track_trueID[itrk] << std::endl;
    } else {
      if (verbosityBASE > 3)
        std::cout << "found rec track: \t" << itrk << "\t track source: " << _track_source[itrk] << " MC: "<< (int)_track_trueID[itrk] << std::endl;
    }
    _mcpart_RecTrackIDs[(int)_track_trueID[itrk]].push_back(itrk);
  }
  Int_t nCurrProj = 0;
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    if (verbosityBASE > 2) std::cout << "current track: " << itrk <<std::endl;
    unsigned short trackSource = _track_source[itrk];
    for(Int_t iproj=nCurrProj; iproj<_nProjections; iproj++){
      if (itrk != _track_ProjTrackID[iproj])
        continue;
      double projectionR = TMath::Sqrt(_track_Proj_x[iproj]*_track_Proj_x[iproj]+_track_Proj_y[iproj]*_track_Proj_y[iproj]);       if(TMath::Abs(_track_Proj_t[iproj])< 2.e-20){
        if (verbosityBASE > 5) std::cout << iproj << "\t projection layer: "<< _track_ProjLayer[iproj] << "\t t: " << _track_Proj_t[iproj] 
                                            << "\t x: " << _track_Proj_x[iproj] << "\t y: " << _track_Proj_y[iproj] << "\t z: " << _track_Proj_z[iproj] << "\t r: " << projectionR << std::endl;

        continue;
      }
      if (verbosityBASE > 2) std::cout << iproj << "\t projection layer: "<< _track_ProjLayer[iproj] << "\t t: " << _track_Proj_t[iproj] 
                                            << "\t x: " << _track_Proj_x[iproj] << "\t y: " << _track_Proj_y[iproj] << "\t z: " << _track_Proj_z[iproj] << "\t r: " << projectionR << std::endl;
      _track_hasTTL[itrk]     = (_track_hasTTL[itrk] || HasTimingLayer(_track_ProjLayer[iproj]));
      if (HasTimingLayer(_track_ProjLayer[iproj])){
        _track_nTTL[itrk]++;
        _track_Proj_Clas[iproj] = 2;      
      }
      if (verbosityBASE > 3) std::cout << "timing layer count: " << _track_nTTL[itrk] << std::endl;
      if (IsTrackerLayer(_track_ProjLayer[iproj])){
        _track_nTrL[itrk]++;
        _track_Proj_Clas[iproj] = 1;
      }
      if (HasFirstTwoLayers(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 4;
      }
      if (IsCaloProjection(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 5;
      }
      if (IsFarForwardProjection(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 6;
      }
      _track_RefProjID[itrk].push_back(iproj);
      nCurrProj = iproj;
    }

    if (trackSource == 1 )
      _track_hasOL[itrk] = (_track_hasOL[itrk] || true);
    if (trackSource == 2 )
      _track_hasIL[itrk]   = (_track_hasIL[itrk] || true);

    if ((int)_track_trueID[itrk] < 0) continue;
    if (_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size() > 1){
      if (verbosityBASE > 2) std::cout << "identified multiple tracks for: " << (int)_track_trueID[itrk] << " current track id: " << itrk << std::endl;
      for (int rtr = 0; rtr < _mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size(); rtr++){
        if (verbosityBASE > 3) std::cout << rtr << "\t MC id " << _mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr) << "\t track source\t " << _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] << std::endl;
        if ( _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] > trackSource &&  _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] == 1 )
          _track_hasOL[itrk] = (_track_hasOL[itrk] || true);
        if ( _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] > trackSource  &&  _track_source[_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].at(rtr)] == 2 )
          _track_hasIL[itrk]   = (_track_hasIL[itrk] || true);
      }
    } else {
      if (verbosityBASE > 2) std::cout << "only one track found for: " << (int)_track_trueID[itrk] << "\t source: " << trackSource<< std::endl;
    }
    if (verbosityBASE > 2) std::cout << "found: " << _track_RefProjID[itrk].size() << "\t projections " << std::endl;
    
  }
}

// **********************************************************************************************
// ****************** cleanup vectors for matching diff rec tracks to MC part ********************
// **********************************************************************************************
void clearMCRecMatchVectors(){
  if (verbosityBASE > 2) std::cout << "clearing MC vectors" << std::endl;
  for(int imc=0; imc<_nMCPart+2 && imc < _maxNMCPart; imc++){
    if (_mcpart_RecTrackIDs[imc].size() > 0){
      _mcpart_RecTrackIDs[imc].clear();
      _mcpart_RecTrackIDs[imc].resize(0);
    }
  }

  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    _track_hasTTL[itrk]   = 0;
    _track_nTTL[itrk]     = 0;
    _track_nTrL[itrk]     = 0;
    _track_hasIL[itrk]    = 0;
    _track_hasOL[itrk]    = 0;
    if (_track_RefProjID[itrk].size() > 0){
      _track_RefProjID[itrk].clear();
      _track_RefProjID[itrk].resize(0);
    }
  }
  for (Int_t iproj=0; iproj<_nProjections; iproj++){
    _track_Proj_Clas[iproj] = 0;
  }
}

