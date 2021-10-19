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
  kBECAL          = 10
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
bool _track_hasTTL[_maxNTracks]          = {0};
bool _track_hasATTL[_maxNTracks]         = {0};
int _track_nTTL[_maxNTracks]             = {0};
bool _track_hasIL[_maxNTracks]           = {0};
bool _track_hasOL[_maxNTracks]           = {0};

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
int* _mcpart_ID                   = new int[_maxNMCPart];
int* _mcpart_ID_parent            = new int[_maxNMCPart];
int* _mcpart_PDG                  = new int[_maxNMCPart];
float* _mcpart_E                  = new float[_maxNMCPart];
float* _mcpart_px                 = new float[_maxNMCPart];
float* _mcpart_py                 = new float[_maxNMCPart];
float* _mcpart_pz                 = new float[_maxNMCPart];
float* _mcpart_Eta                = new float[_maxNMCPart];
float* _mcpart_Phi                = new float[_maxNMCPart];
std::vector<int> _mcpart_RecTrackIDs[_maxNMCPart];


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
const int _maxProjectionLayers    = 8;
Int_t layerIndexForward[9]  =  {20, 21, 22, 23, 24, 0, 1, 2, -1};
void ResetLayerIndexForward(){
  if (!_is_ALLSILICON){
    layerIndexForward[0]  = 50;
    layerIndexForward[1]  = 51;
    layerIndexForward[2]  = 52;
    layerIndexForward[3]  = 53;
    layerIndexForward[4]  = 54;
  }
}
TString GetProjectionNameFromIndex(int projindex)
{
  switch (projindex){
    case 0:     return "FTTL_0";
    case 1:     return "FTTL_1";
    case 2:     return "FTTL_2";
    case 3:     return "ETTL_0";
    case 4:     return "ETTL_1";
    case 5:    return "FHCAL";
    case 6:    return "FEMC";
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
            
    case 40:    return "CENTAL_0"; //"BARREL_0";
    case 41:    return "CENTAL_1"; //"BARREL_1";
    case 42:    return "CENTAL_2"; //"BARREL_2";
    case 43:    return "CENTAL_3"; //"BARREL_3";
    case 44:    return "CENTAL_4"; //"BARREL_4";
    case 45:    return "CENTAL_5"; //"BARREL_5";

    case 50:    return "FORWARD_0";
    case 51:    return "FORWARD_1";
    case 52:    return "FORWARD_2";
    case 53:    return "FORWARD_3";
    case 54:    return "FORWARD_4";
    case 55:    return "FORWARD_5";
  
    default:   return "NOTHING";
  }
}

Int_t ReturnIndexForwardLayer(Int_t layerID){
  if (_is_ALLSILICON){
    switch (layerID){
      case 20: return 0;
      case 21: return 1;
      case 22: return 2;
      case 23: return 3;
      case 24: return 4;
      case 0: return 5;
      case 1: return 6;
      case 2: return 7;
      default: return -1;  
    }
  } else {
    switch (layerID){
      case 50: return 0;
      case 51: return 1;
      case 52: return 2;
      case 53: return 3;
      case 54: return 4;
      case 0: return 5;
      case 1: return 6;
      case 2: return 7;
      default: return -1;  
    }
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
      case 50: 
      case 51: 
      case 40:
      case 41:
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
      case 50: 
      case 51: 
      case 52: 
      case 53: 
      case 54: 
      case 55:   
      case 40:
      case 41:
      case 42:
      case 43:
      case 44:
      case 45:
        return kTRUE;
      default:
        return kFALSE;
    }    
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

Bool_t HasTimingLayerAfterECal(Int_t layerID){
  switch (layerID){
    case 2:
    case 8:
      return kTRUE;
    default:
      return kFALSE;
  }
  return kFALSE;
}

int GetCorrectMCArrayEntry(float objectTrueID){
  for(Int_t imc=0; imc<_nMCPart; imc++){
    if(objectTrueID==_mcpart_ID[imc]){
      return imc;
    }
  }
  return -1;
}
