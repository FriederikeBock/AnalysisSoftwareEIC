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
const int _maxNProjections = 20000;
const int _maxNMCPart = 100000;
const int _maxNHepmcp = 100000;
int verbosityBASE = 0;
bool _useAlternateForProjections = false;
float _nEventsTree;
// TOF resolution
 double sigmat             = 25e-3; // ns default ECCE
// double sigmat             = 30e-3; // ns
// double sigmat             = 35e-3; // ns
// double sigmat             = 40e-3; // ns
// double sigmat             = 50e-3; // ns
// random number generator
TRandom3 r3;
bool isZeroField      = false;

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

struct projStrct{
  projStrct(): caloid(-1), eta_Calo(-1000), phi_Calo(-1000) {}
  int caloid;
  float eta_Calo;
  float phi_Calo;
} ;


struct matchingCalStrct{
  matchingCalStrct(): caloid(-1), id(-1), dPhi(-1000), dEta(-1000), dX(-1000), dY(-1000) {}
  int caloid;
  int id;
  float dPhi;
  float dEta;
  float dX;
  float dY;
} ;

struct tofStrct{
  tofStrct(): pathlength(0), recTime(0), genTime(0) {}
  tofStrct(float pl, float rT, float gT, bool iPE){
    pathlength  = pl;
    recTime     = rT;
    genTime     = gT;
    isProbEScat = iPE;
  }
  float pathlength;
  float recTime;
  float genTime;
  bool isProbEScat;
} ;



bool caloEnabled[20]      = {0};
bool tracksEnabled        = 0;
bool vertexEnabled        = 0;
bool extendTimingInfo     = 0;
bool hepMCavailable       = 0;
bool xSectionEnabled      = 0;

// Event level info
float _cross_section;
float _event_weight;
int _n_generator_accepted;

// track hits
int _nHitsLayers;
int* _hits_layerID         = new int[_maxNHits];
float* _hits_x             = new float[_maxNHits];
float* _hits_y             = new float[_maxNHits];
float* _hits_z             = new float[_maxNHits];
float* _hits_t             = new float[_maxNHits];
float* _hits_edep          = new float[_maxNHits];

// towers
int _nTowers_CEMC;
float* _tower_CEMC_E          = new float[_maxNTowersCentral];
int* _tower_CEMC_iEta         = new int[_maxNTowersCentral];
int* _tower_CEMC_iPhi         = new int[_maxNTowersCentral];
int* _tower_CEMC_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EEMC;
float* _tower_EEMC_E          = new float[_maxNTowersCentral];
int* _tower_EEMC_iEta         = new int[_maxNTowersCentral];
int* _tower_EEMC_iPhi         = new int[_maxNTowersCentral];
int* _tower_EEMC_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EEMCG;
float* _tower_EEMCG_E          = new float[_maxNTowersCentral];
int* _tower_EEMCG_iEta         = new int[_maxNTowersCentral];
int* _tower_EEMCG_iPhi         = new int[_maxNTowersCentral];
int* _tower_EEMCG_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_EHCAL;
float* _tower_EHCAL_E          = new float[_maxNTowersCentral];
int* _tower_EHCAL_iEta         = new int[_maxNTowersCentral];
int* _tower_EHCAL_iPhi         = new int[_maxNTowersCentral];
int* _tower_EHCAL_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALIN;
float* _tower_HCALIN_E          = new float[_maxNTowersCentral];
int* _tower_HCALIN_iEta         = new int[_maxNTowersCentral];
int* _tower_HCALIN_iPhi         = new int[_maxNTowersCentral];
int* _tower_HCALIN_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_HCALOUT;
float* _tower_HCALOUT_E          = new float[_maxNTowersCentral];
int* _tower_HCALOUT_iEta         = new int[_maxNTowersCentral];
int* _tower_HCALOUT_iPhi         = new int[_maxNTowersCentral];
int* _tower_HCALOUT_trueID       = new int[_maxNTowersCentral];

// towers
int _nTowers_DRCALO;
float* _tower_DRCALO_E          = new float[_maxNTowersDR];
int* _tower_DRCALO_iEta         = new int[_maxNTowersDR];
int* _tower_DRCALO_iPhi         = new int[_maxNTowersDR];
int* _tower_DRCALO_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_FOCAL;
float* _tower_FOCAL_E          = new float[_maxNTowersDR];
int* _tower_FOCAL_iEta         = new int[_maxNTowersDR];
int* _tower_FOCAL_iPhi         = new int[_maxNTowersDR];
int* _tower_FOCAL_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_LFHCAL;
float* _tower_LFHCAL_E          = new float[_maxNTowersDR];
int* _tower_LFHCAL_iEta         = new int[_maxNTowersDR];
int* _tower_LFHCAL_iPhi         = new int[_maxNTowersDR];
int* _tower_LFHCAL_iL           = new int[_maxNTowersDR];
int* _tower_LFHCAL_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_BECAL;
float* _tower_BECAL_E          = new float[_maxNTowersDR];
int* _tower_BECAL_iEta         = new int[_maxNTowersDR];
int* _tower_BECAL_iPhi         = new int[_maxNTowersDR];
int* _tower_BECAL_trueID       = new int[_maxNTowersDR];

// towers
int _nTowers_FHCAL;
float* _tower_FHCAL_E          = new float[_maxNTowers];
int* _tower_FHCAL_iEta         = new int[_maxNTowers];
int* _tower_FHCAL_iPhi         = new int[_maxNTowers];
int* _tower_FHCAL_trueID       = new int[_maxNTowers];

// towers
int _nTowers_FEMC;
float* _tower_FEMC_E           = new float[_maxNTowers];
int* _tower_FEMC_iEta          = new int[_maxNTowers];
int* _tower_FEMC_iPhi          = new int[_maxNTowers];
int* _tower_FEMC_trueID        = new int[_maxNTowers];

// clusters
int _nclusters_FHCAL;
float* _clusters_FHCAL_E          = new float[_maxNclusters];
float* _clusters_FHCAL_Eta        = new float[_maxNclusters];
float* _clusters_FHCAL_Phi        = new float[_maxNclusters];
int* _clusters_FHCAL_NTower       = new int[_maxNclusters];
int* _clusters_FHCAL_trueID       = new int[_maxNclusters];

// clusters
int _nclusters_FEMC;
float* _clusters_FEMC_E           = new float[_maxNclusters];
float* _clusters_FEMC_Eta         = new float[_maxNclusters];
float* _clusters_FEMC_Phi         = new float[_maxNclusters];
int* _clusters_FEMC_NTower        = new int[_maxNclusters];
int* _clusters_FEMC_trueID        = new int[_maxNclusters];

// clusters
int _nclusters_HCALOUT;
float* _clusters_HCALOUT_E           = new float[_maxNclusters];
float* _clusters_HCALOUT_Eta         = new float[_maxNclusters];
float* _clusters_HCALOUT_Phi         = new float[_maxNclusters];
int* _clusters_HCALOUT_NTower        = new int[_maxNclusters];
int* _clusters_HCALOUT_trueID        = new int[_maxNclusters];

// clusters
int _nclusters_HCALIN;
float* _clusters_HCALIN_E           = new float[_maxNclusters];
float* _clusters_HCALIN_Eta         = new float[_maxNclusters];
float* _clusters_HCALIN_Phi         = new float[_maxNclusters];
int* _clusters_HCALIN_NTower        = new int[_maxNclusters];
int* _clusters_HCALIN_trueID        = new int[_maxNclusters];

// vertex
float _vertex_x;
float _vertex_y;
float _vertex_z;
float _vertex_true_x;
float _vertex_true_y;
float _vertex_true_z;

Int_t trackID_ScatteredElectron = -1;

// tracks
int _nTracks;
float* _track_ID                 = new float[_maxNTracks];
float* _track_trueID             = new float[_maxNTracks];
float* _track_px                 = new float[_maxNTracks];
float* _track_py                 = new float[_maxNTracks];
float* _track_pz                 = new float[_maxNTracks];
float* _track_x                  = new float[_maxNTracks];
float* _track_y                  = new float[_maxNTracks];
float* _track_z                  = new float[_maxNTracks];
short* _track_charge             = new short[_maxNTracks];
float* _track_dca2D              = new float[_maxNTracks];
float* _track_dca                = new float[_maxNTracks];
float* _track_chi2               = new float[_maxNTracks];
float* _track_ndf                = new float[_maxNTracks];
unsigned short* _track_source             = new unsigned short[_maxNTracks];

std::array<std::vector<int>, _maxNTracks> _track_RefProjID;
// Initializes all elements to 0
std::array<bool, _maxNTracks> _track_hasTTL{{}};
std::array<int, _maxNTracks> _track_nTTL{{}};
std::array<int, _maxNTracks> _track_nTrL{{}};
std::array<bool, _maxNTracks> _track_hasIL{{}};
std::array<bool, _maxNTracks> _track_hasOL{{}};
std::array<matchingCalStrct, _maxNTracks> _track_matchECal{{}};
std::array<matchingCalStrct, _maxNTracks> _track_matchHCal{{}};
std::array<std::vector<tofStrct>, _maxNTracks> _track_TOFmeas;

// hadron PID
float* _track_pion_LL           = new float[_maxNTracks];
float* _track_kaon_LL           = new float[_maxNTracks];
float* _track_proton_LL         = new float[_maxNTracks];

int _nProjections;
float* _track_ProjTrackID                 = new float[_maxNProjections];
int* _track_ProjLayer                 = new int[_maxNProjections];
float* _track_Proj_x             = new float[_maxNProjections];
float* _track_Proj_y             = new float[_maxNProjections];
float* _track_Proj_z             = new float[_maxNProjections];
float* _track_Proj_t             = new float[_maxNProjections];
float* _track_Proj_px             = new float[_maxNProjections];
float* _track_Proj_py             = new float[_maxNProjections];
float* _track_Proj_pz             = new float[_maxNProjections];
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
float* _mcpart_x                 = new float[_maxNMCPart];
float* _mcpart_y                 = new float[_maxNMCPart];
float* _mcpart_z                 = new float[_maxNMCPart];
float* _mcpart_Eta                = new float[_maxNMCPart];
float* _mcpart_Phi                = new float[_maxNMCPart];
int* _mcpart_BCID                = new int[_maxNMCPart];
std::array<std::vector<int>, _maxNMCPart> _mcpart_RecTrackIDs;
std::array<std::vector<projStrct>, _maxNMCPart> _mcpart_EcalProjs;
std::array<std::vector<projStrct>, _maxNMCPart> _mcpart_HcalProjs;

int _nHepmcp;
int _hepmcp_procid;
float _hepmcp_vtx_x;
float _hepmcp_vtx_y;
float _hepmcp_vtx_z;
float _hepmcp_vtx_t;
float _hepmcp_x1;
float _hepmcp_x2;
float _hepmcp_Q2;
int*  _hepmcp_BCID            = new int[_maxNHepmcp];
// float*  _hepmcp_ID_parent = new float[_maxNHepmcp];
int*  _hepmcp_status          = new int[_maxNHepmcp];
int*  _hepmcp_PDG             = new int[_maxNHepmcp];
float*  _hepmcp_E             = new float[_maxNHepmcp];
float*  _hepmcp_px            = new float[_maxNHepmcp];
float*  _hepmcp_py            = new float[_maxNHepmcp];
float*  _hepmcp_pz            = new float[_maxNHepmcp];
int*  _hepmcp_m1              = new int[_maxNHepmcp];
int*  _hepmcp_m2              = new int[_maxNHepmcp];

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
      inputTree->SetBranchAddress("tracks_x",             _track_x);
      inputTree->SetBranchAddress("tracks_y",             _track_y);
      inputTree->SetBranchAddress("tracks_z",             _track_z);
      inputTree->SetBranchAddress("tracks_trueID",        _track_trueID);
      inputTree->SetBranchAddress("tracks_source",        _track_source);
      inputTree->SetBranchAddress("tracks_charge",        _track_charge);
      inputTree->SetBranchAddress("tracks_dca",           _track_dca);
      inputTree->SetBranchAddress("tracks_dca_2d",        _track_dca2D);
      inputTree->SetBranchAddress("tracks_ndf",           _track_ndf);
      inputTree->SetBranchAddress("tracks_chi2",          _track_chi2);
      
      inputTree->SetBranchAddress("track_pion_LL",	  _track_pion_LL);
      inputTree->SetBranchAddress("track_kaon_LL",	  _track_kaon_LL);
      inputTree->SetBranchAddress("track_proton_LL",	  _track_proton_LL);
 
      inputTree->SetBranchAddress("nProjections",         &_nProjections);
      inputTree->SetBranchAddress("track_ProjTrackID",    _track_ProjTrackID);
      inputTree->SetBranchAddress("track_ProjLayer",      _track_ProjLayer);

      inputTree->SetBranchAddress("track_TLP_x",           _track_Proj_x);
      inputTree->SetBranchAddress("track_TLP_y",           _track_Proj_y);
      inputTree->SetBranchAddress("track_TLP_z",           _track_Proj_z);
      inputTree->SetBranchAddress("track_TLP_t",           _track_Proj_t);
      inputTree->SetBranchAddress("track_TLP_px",          _track_Proj_px);
      inputTree->SetBranchAddress("track_TLP_py",          _track_Proj_py);
      inputTree->SetBranchAddress("track_TLP_pz",          _track_Proj_pz);
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
    if (caloEnabled[kHCALOUT]){
      if (inputTree->GetBranchStatus("cluster_HCALOUT_N") ){
        inputTree->SetBranchAddress("cluster_HCALOUT_N",                 &_nclusters_FEMC);
        inputTree->SetBranchAddress("cluster_HCALOUT_E",                 _clusters_FEMC_E);
        inputTree->SetBranchAddress("cluster_HCALOUT_Eta",              _clusters_FEMC_Eta);
        inputTree->SetBranchAddress("cluster_HCALOUT_Phi",              _clusters_FEMC_Phi);
        inputTree->SetBranchAddress("cluster_HCALOUT_NTower",              _clusters_FEMC_NTower);
        inputTree->SetBranchAddress("cluster_HCALOUT_trueID",            _clusters_FEMC_trueID);
      }
    }
    if (caloEnabled[kHCALIN]){
      if (inputTree->GetBranchStatus("cluster_HCALIN_N") ){
        inputTree->SetBranchAddress("cluster_HCALIN_N",                 &_nclusters_FEMC);
        inputTree->SetBranchAddress("cluster_HCALIN_E",                 _clusters_FEMC_E);
        inputTree->SetBranchAddress("cluster_HCALIN_Eta",              _clusters_FEMC_Eta);
        inputTree->SetBranchAddress("cluster_HCALIN_Phi",              _clusters_FEMC_Phi);
        inputTree->SetBranchAddress("cluster_HCALIN_NTower",            _clusters_FEMC_NTower);
        inputTree->SetBranchAddress("cluster_HCALIN_trueID",            _clusters_FEMC_trueID);
      }
    }
    // vertex

    if (inputTree->GetBranchStatus("vertex_x") ){
      vertexEnabled = 1;
      inputTree->SetBranchAddress("vertex_x",                     &_vertex_x);
      inputTree->SetBranchAddress("vertex_y",                     &_vertex_y);
      inputTree->SetBranchAddress("vertex_z",                     &_vertex_z);
      inputTree->SetBranchAddress("vertex_true_x",                &_vertex_true_x);
      inputTree->SetBranchAddress("vertex_true_y",                &_vertex_true_y);
      inputTree->SetBranchAddress("vertex_true_z",                &_vertex_true_z);
    }
    // MC particles

    inputTree->SetBranchAddress("nMCPart",          &_nMCPart);
    inputTree->SetBranchAddress("mcpart_ID",        _mcpart_ID);
    inputTree->SetBranchAddress("mcpart_ID_parent", _mcpart_ID_parent);
    inputTree->SetBranchAddress("mcpart_PDG",       _mcpart_PDG);
    inputTree->SetBranchAddress("mcpart_E",         _mcpart_E);
    inputTree->SetBranchAddress("mcpart_px",        _mcpart_px);
    inputTree->SetBranchAddress("mcpart_py",        _mcpart_py);
    inputTree->SetBranchAddress("mcpart_pz",        _mcpart_pz);
    inputTree->SetBranchAddress("mcpart_BCID",      _mcpart_BCID);
    if (inputTree->GetBranchStatus("mcpart_x") ){
      extendTimingInfo = 1;
      inputTree->SetBranchAddress("mcpart_x",         _mcpart_x);
      inputTree->SetBranchAddress("mcpart_y",         _mcpart_y);
      inputTree->SetBranchAddress("mcpart_z",         _mcpart_z);
    }
    

    if (inputTree->GetBranchStatus("nHepmcp") ){
      hepMCavailable = 1;
      inputTree->SetBranchAddress("nHepmcp", &_nHepmcp);
      inputTree->SetBranchAddress("hepmcp_procid", &_hepmcp_procid);
      inputTree->SetBranchAddress("hepmcp_vtx_x", &_hepmcp_vtx_x);
      inputTree->SetBranchAddress("hepmcp_vtx_y", &_hepmcp_vtx_y);
      inputTree->SetBranchAddress("hepmcp_vtx_z", &_hepmcp_vtx_z);
      inputTree->SetBranchAddress("hepmcp_vtx_t", &_hepmcp_vtx_t);
      inputTree->SetBranchAddress("hepmcp_x1", &_hepmcp_x1);
      inputTree->SetBranchAddress("hepmcp_x2", &_hepmcp_x2);
      inputTree->SetBranchAddress("hepmcp_Q2", &_hepmcp_Q2);

      //    inputTree->SetBranchAddress("hepmcp_ID_parent", _hepmcp_ID_parent, "hepmcp_ID_parent[nHepmcp]/F");
      inputTree->SetBranchAddress("hepmcp_status", _hepmcp_status);
      inputTree->SetBranchAddress("hepmcp_PDG", _hepmcp_PDG);
      inputTree->SetBranchAddress("hepmcp_E", _hepmcp_E);
      inputTree->SetBranchAddress("hepmcp_px", _hepmcp_px);
      inputTree->SetBranchAddress("hepmcp_py", _hepmcp_py);
      inputTree->SetBranchAddress("hepmcp_pz", _hepmcp_pz);
      inputTree->SetBranchAddress("hepmcp_BCID", _hepmcp_BCID);
      inputTree->SetBranchAddress("hepmcp_m1", _hepmcp_m1);
      inputTree->SetBranchAddress("hepmcp_m2", _hepmcp_m2);
    }
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
const int _maxProjectionLayers    = 36;
Int_t layerIndexHist[36]  =  { 0, 1, 2, 3, 4, 7, 8,        // TTLs
                                 50, 51, 52, 53, 54,      // FST
                                 113, 114, 115, 116,      // EFST
                                 155, 156, 157, 150, 151, // SVTX & BARR
                                 110, 111, 112, 120, 130,  // mu RWell 
                                 6, 60, 61, 62, 63, 66, 67, // calo's 
                                 35, 36, 37 // cherenkov PID
};
Int_t enableLayerHist[36]  =  { 1, 0, 0, 1, 0, 1, 0,        // TTLs
                                 1, 1, 1, 1, 1,      // FST
                                 1, 1, 1, 1,      // EFST
                                 1, 1, 1, 1, 1, // SVTX & BARR
                                 1, 1, 1, 0, 0,  // mu RWell 
                                 1, 1, 1, 1, 1, 1, 1, // calo's 
                                 1, 1, 1 // cherenkov PID
};

//__________________________________________________________________________________________________________
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
    // EST
    case 113:    return "EST_0";
    case 114:    return "EST_1";
    case 115:    return "EST_2";
    case 116:    return "EST_3";
    case 117:    return "EST_4";
    case 118:    return "EST_5";
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
    case 6:    return "FEMC";
    case 60:   return "EHCAL";
    case 61:   return "EEMC";
    case 62:   return "HCALIN";
    case 63:   return "HCALOUT";
    case 64:   return "CEMC";
    case 65:   return "EEMC_glass_0";
    case 66:   return "BECAL";
    case 5:    //return "FHCAL";
    case 67:   
      return "LFHCAL";
    case 140:  return "LFHCAL_0";
    case 141:  return "LFHCAL_1";
    case 142:  return "LFHCAL_2";
    case 143:  return "LFHCAL_3";
    case 144:  return "LFHCAL_4";
    case 145:  return "LFHCAL_5";
    case 146:  return "LFHCAL_6";
    case 147:  return "LFHCAL_7";
    
    case 35: return "hpDIRC";
    case 36: return "mRICH";
    case 37: return "dRICH";
    
    default:   return "NOTHING";
  }
}

//__________________________________________________________________________________________________________
Int_t GetRegionFromEta(float etaMC){
  if (etaMC > nominalEtaRegion[0][0] && etaMC < nominalEtaRegion[0][1] )
    return 0;
  else if (etaMC > nominalEtaRegion[1][0] && etaMC < nominalEtaRegion[1][1] )
    return 1;
  else if (etaMC > nominalEtaRegion[2][0] && etaMC < nominalEtaRegion[2][1] )
    return 2;
  else 
    return -1;
}

//__________________________________________________________________________________________________________
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
    
    // cherenkov PID
    case 35: return 1;    // hpDIRC
    case 36: return 0;    // mRICH
    case 37: return 2;    // dRICH
    
    default:   return 0;
  }
}

//__________________________________________________________________________________________________________
int ReturnCalorimeterFromProjectionIndex(int projID){
  switch (projID){
    case 1: return kDRCALO;
//     case 5: return kFHCAL;
    case 6: return kFEMC;
    case 60: return kEHCAL;
    case 61: return kEEMC;
    case 62: return kHCALIN;
    case 63: return kHCALOUT;
    case 64: return kCEMC;
    case 65: return kEEMCG;
    case 67: 
    case 5: 
      return kLFHCAL;
    case 66: return kBECAL;
    case 85: return kFOCAL;
    default:
      // std::cout << "ReturnCalorimeterFromProjectionIndex: projID " << projID << " not defined, returning -1" << std::endl;
      return -1;
  }
  return -1;
}

//__________________________________________________________________________________________________________
Int_t ReturnIndexForwardLayer(Int_t layerID){
  for (Int_t i = 0; i < _maxProjectionLayers; i++){
    if (layerIndexHist[i] == layerID)
      return i;
  }
  return -1;
}

//__________________________________________________________________________________________________________
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
      case 113: // EST_0
      case 114: // EST_1
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

//__________________________________________________________________________________________________________
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
      // EFST
      case 113:
      case 114:
      case 115:
      case 116:
      case 117:
      case 118:
      case 119:
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

//__________________________________________________________________________________________________________
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

//__________________________________________________________________________________________________________
Bool_t IsECalProjection(Int_t layerID, bool alternate = false){
  if (!alternate){
    switch (layerID){  
      case 6: 
      case 61: 
      case 64:
      case 65:
      case 66:
        return kTRUE;
      default:
        return kFALSE;
    }
  } else {
    switch (layerID){      
      case 0:     
      case 3:     
      case 7:     
        return kTRUE;
      default:
        return kFALSE;
    }
  }
  return kFALSE;
}

//__________________________________________________________________________________________________________
Bool_t IsHCalProjection(Int_t layerID, bool alternate = false){
  if (!alternate){
    switch (layerID){  
      case 5: 
      case 60: 
      case 62:
      case 63:
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
  } else {
    switch (layerID){      
      case 0:     
      case 3:     
      case 7:     
        return kTRUE;
      default:
        return kFALSE;
    }
  }
  return kFALSE;
}

//__________________________________________________________________________________________________________
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

//__________________________________________________________________________________________________________
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

//__________________________________________________________________________________________________________
bool  GetMaxCoordinateCalo( float &maxcoord, int caloID){
  switch (caloID){
    case kDRCALO: 
      maxcoord = 262.+ 15.;
      return kFALSE;
    case kFHCAL: 
      maxcoord = 262.+ 15.;
      return kFALSE;
    case kFEMC: 
      maxcoord = 182.655 + 15.;
      return kFALSE;
    case kEHCAL: 
      maxcoord = 260 + 15;
      return kFALSE;
    case kEEMC: 
      maxcoord = 61. + 6.;
      return kFALSE;
    case kHCALIN:
      maxcoord = 230;
      return kTRUE;
    case kHCALOUT:
      maxcoord = 350;
      return kTRUE;
    case kCEMC: 
      maxcoord = 230;
      return kTRUE;
    case kEEMCG: 
      maxcoord = 67.;
      return kFALSE;
    case kLFHCAL: 
      maxcoord = 262.+ 15.;
      return kFALSE;
    case kBECAL: 
      maxcoord = 230;
      return kTRUE;
    case kFOCAL:
      maxcoord = 262.+ 15.;
      return kFALSE;
    default:
      return kFALSE;
  }
  return kFALSE;
}

//__________________________________________________________________________________________________________
bool  GetMinCoordinateCalo( float &mincoord, int caloID){
  switch (caloID){
    case kDRCALO: 
      mincoord = 262.+ 15.;
      return kFALSE;
    case kFHCAL: 
      mincoord = 262.+ 15.;
      return kFALSE;
    case kFEMC: 
      mincoord = 182.655 + 15.;
      return kFALSE;
    case kEHCAL: 
      mincoord = 260 + 15;
      return kFALSE;
    case kEEMC: 
      mincoord = 61. + 6.;
      return kFALSE;
    case kHCALIN:
      mincoord = -305;
      return kTRUE;
    case kHCALOUT:
      mincoord = -305;
      return kTRUE;
    case kCEMC: 
      mincoord = -305;
      return kTRUE;
    case kEEMCG: 
      mincoord = 67.;
      return kFALSE;
    case kLFHCAL: 
      mincoord = 262.+ 15.;
      return kFALSE;
    case kBECAL: 
      mincoord = -305;
      return kTRUE;
    case kFOCAL:
      mincoord = 262.+ 15.;
      return kFALSE;
    default:
      return kFALSE;
  }
  return kFALSE;
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
  trackID_ScatteredElectron = -1;
  for(Int_t itrk=0; itrk<_nTracks; itrk++){
    if (verbosityBASE > 2) std::cout << "current track: " << itrk <<std::endl;
    if (verbosityBASE > 5){
      TVector3 recpartvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
      TVector3 mcpartvec(_mcpart_px[(int)_track_trueID[itrk]],_mcpart_py[(int)_track_trueID[itrk]],_mcpart_pz[(int)_track_trueID[itrk]]);
      float receta = recpartvec.Eta();
      
      float trueeta = mcpartvec.Eta();
      std::cout << "\tTRKTRUEID " << (int)_track_trueID[itrk] << "\tPDG: " << _mcpart_PDG[(int)_track_trueID[itrk]] << "\tpT: " <<   TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)) << "\tEta " << receta << "\t track source\t" << _track_source[itrk] << " MC: "<< (int)_track_trueID[itrk] << " ndf: " << _track_ndf[itrk] << " chi2: " << _track_chi2[itrk]<< std::endl;
    }
    unsigned short trackSource = _track_source[itrk];
    TVector3 trackVecP(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
    if (_track_trueID[itrk] >= 0 && trackSource == 0){
      Int_t mcID = _track_trueID[itrk];
      if (_mcpart_PDG[mcID] == 11){
        int bcid = _mcpart_BCID[mcID]-1;
        if (_hepmcp_status[bcid] == 1 && _hepmcp_status[_hepmcp_m1[bcid]-1] > 1)
          trackID_ScatteredElectron = itrk;
      }      
    }
    
    for(Int_t iproj=nCurrProj; iproj<_nProjections; iproj++){
      if (itrk != _track_ProjTrackID[iproj])
        continue;
//       if (_track_Proj_t[iproj] < 0 )
// //         continue;
      double projectionR = TMath::Sqrt(_track_Proj_x[iproj]*_track_Proj_x[iproj]+_track_Proj_y[iproj]*_track_Proj_y[iproj]);      
      double hit3d = TMath::Sqrt(_track_Proj_true_x[iproj]*_track_Proj_true_x[iproj]+_track_Proj_true_y[iproj]*_track_Proj_true_y[iproj]+_track_Proj_true_z[iproj]*_track_Proj_true_z[iproj]);      
      if(TMath::Abs(_track_Proj_true_t[iproj]) == 0. || _track_Proj_true_t[iproj] == -10000.){
//         if (verbosityBASE > 5) std::cout << "failed: " << iproj << "\t projection layer: "<< _track_ProjLayer[iproj] << "\t t: " << _track_Proj_t[iproj] 
//                                          << "\t x: " << _track_Proj_x[iproj] << "\t y: " << _track_Proj_y[iproj] << "\t z: " << _track_Proj_z[iproj] << "\t r: " 
//                                          << projectionR << "\t true x " << _track_Proj_true_x[iproj] << "\t true y " << _track_Proj_true_y[iproj] << "\t true z " << _track_Proj_true_z[iproj]<< "\t true t " << _track_Proj_true_t[iproj] <<  std::endl;

        continue;
      }
      if (verbosityBASE > 2) std::cout << iproj << "\t projection layer: "<< _track_ProjLayer[iproj] << "\t t: " << _track_Proj_t[iproj] 
                                       << "\t x: " << _track_Proj_x[iproj] << "\t y: " << _track_Proj_y[iproj] << "\t z: " << _track_Proj_z[iproj] << "\t r: " 
                                       << projectionR << "\t true x " << _track_Proj_true_x[iproj] << "\t true y " << _track_Proj_true_y[iproj] << "\t true z " << _track_Proj_true_z[iproj]<< "\t true t " << _track_Proj_true_t[iproj] <<  std::endl;
      if (HasTimingLayer(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 2;      
      }
      if (IsTrackerLayer(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 1;
      }
      if (HasFirstTwoLayers(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 4;
      }
      if (IsCaloProjection(_track_ProjLayer[iproj])){
        float minVal = 0;
        float maxVal = 0;
        int calID    = ReturnCalorimeterFromProjectionIndex(_track_ProjLayer[iproj]); 
        bool useZ = GetMinCoordinateCalo( minVal, calID);
        GetMaxCoordinateCalo( maxVal, calID);
        if (useZ){
          if ( !(_track_Proj_z[iproj] > minVal && _track_Proj_z[iproj] < maxVal)) continue;
        } else {
          
          if ( projectionR > maxVal ) continue;
        } 
        _track_Proj_Clas[iproj] = 5;
      }
      if (IsFarForwardProjection(_track_ProjLayer[iproj])){
        _track_Proj_Clas[iproj] = 6;
      }
      if (verbosityBASE > 2) std::cout << iproj << "\t projection layer: "<< _track_ProjLayer[iproj] << "\t class: " << _track_Proj_Clas[iproj] << std::endl;
      if (_track_Proj_Clas[iproj] < 5){
        if (IsTrackerLayer(_track_ProjLayer[iproj]) && _track_Proj_true_t[iproj] != 0 && hit3d > 0.1){
          
          _track_nTrL[itrk]++;
        }
        if (HasTimingLayer(_track_ProjLayer[iproj]) && _track_Proj_true_t[iproj] != 0 && hit3d > 0.1){
          
          // correct path length for geometry and vertex effects
          float corrPathTTL       = TMath::Sqrt(pow(_track_Proj_true_x[iproj],2) + pow(_track_Proj_true_y[iproj],2) + pow(_track_Proj_true_z[iproj],2)) - 
                                    TMath::Sqrt(pow(_track_Proj_x[iproj],2) + pow( _track_Proj_y[iproj],2) + pow(_track_Proj_z[iproj],2));
          float corrPathVtx       = 0;
          if ((int)_track_trueID[itrk] > -1 && extendTimingInfo){
            // point of closest approach - MC origin
            corrPathVtx            = TMath::Sqrt(pow(_track_x[itrk],2) + pow( _track_y[itrk],2) + pow(_track_z[itrk],2)) - 
                                     TMath::Sqrt(pow(_mcpart_x[(int)_track_trueID[itrk]],2) + pow(_mcpart_y[(int)_track_trueID[itrk]],2) + pow(_mcpart_z[(int)_track_trueID[itrk]],2));
          }
          Float_t pathlengthcorr  = _track_Proj_t[iproj] + corrPathTTL + corrPathVtx; // 
          
          if (TMath::Abs(corrPathTTL) > 5 && verbosityBASE > 3){
            float_t Rproj = TMath::Sqrt(_track_Proj_x[iproj]*_track_Proj_x[iproj] + _track_Proj_y[iproj]*_track_Proj_y[iproj]  ) ;
            float_t Rtrue = TMath::Sqrt(_track_Proj_true_x[iproj]*_track_Proj_true_x[iproj] + _track_Proj_true_y[iproj]*_track_Proj_true_y[iproj] ) ;
            std::cout << "layer ID: "<< _track_ProjLayer[iproj] << "\t p\t" << trackVecP.Mag() << "\t path length org: " << _track_Proj_t[iproj] << "\t TTL corr \t" <<  corrPathTTL << "\t Vtx corr \t"<< corrPathVtx << " corrected path length "<< pathlengthcorr << " \t eta \t " << trackVecP.Eta() << std::endl;
          
            std::cout <<  "\t\t projection: " << _track_Proj_x[iproj] << "\t" << _track_Proj_y[iproj] << "\t" <<  _track_Proj_z[iproj] << "\t R " << Rproj << 
                          "\t\t hit: " << _track_Proj_true_x[iproj] << "\t" << _track_Proj_true_y[iproj] << "\t" <<  _track_Proj_true_z[iproj]  << "\t" <<  _track_Proj_true_t[iproj] << "\t R \t" << Rtrue << "\n\n"<< std::endl;
          }
          if ( TMath::Abs(corrPathVtx) > 3 && verbosityBASE > 3){
            float_t Rproj = TMath::Sqrt(pow(_track_x[itrk],2) + pow( _track_y[itrk],2) + pow(_track_z[itrk],2));
            float_t Rtrue = 0;
            if ((int)_track_trueID[itrk] > -1 && extendTimingInfo){
              Rtrue = TMath::Sqrt(pow(_mcpart_x[(int)_track_trueID[itrk]],2) + pow(_mcpart_y[(int)_track_trueID[itrk]],2) + pow(_mcpart_z[(int)_track_trueID[itrk]],2));
              std::cout << "layer ID: "<< _track_ProjLayer[iproj] << "\t p\t" << trackVecP.Mag() << "\t path length org: " << _track_Proj_t[iproj] << "\t TTL corr \t" <<  corrPathTTL << "\t Vtx corr \t"<< corrPathVtx << " corrected path length "<< pathlengthcorr << " \t eta \t " << trackVecP.Eta() << "\n\n" << std::endl;
              std::cout <<  "\t\t track start: " << _track_x[itrk] << "\t" << _track_y[itrk] << "\t" <<  _track_z[itrk] << 
                          "\t\t true start: " << _mcpart_x[(int)_track_trueID[itrk]] << "\t" << _mcpart_y[(int)_track_trueID[itrk]] << "\t" <<  _mcpart_z[(int)_track_trueID[itrk]]  << "\t 3D R\t" << Rproj << "\t" << Rtrue << "\n\n"<< std::endl;
            }
          }
          // check whether this is an electron in the backward or barrel with momentum > 3 GeV => assume that's the scattered electron
          bool isProbEScat = false;
          if ((int)_track_trueID[itrk] > -1){
            if ( (fabs(_mcpart_PDG[(int)_track_trueID[itrk]]) == 11 && trackVecP.Eta() < 0.5 && trackVecP.Mag() > 1))
              isProbEScat = true;
          }
          // only store tof info if vertex correction is smaller than 5cm, otherwise this is a secondary or something went terribly wrong otherwise
          if (TMath::Abs(corrPathVtx) < 5){ 
            _track_nTTL[itrk]++;
            _track_hasTTL[itrk]     = (_track_hasTTL[itrk] || HasTimingLayer(_track_ProjLayer[iproj]));

            // tof construct
            tofStrct tempTOF = tofStrct(pathlengthcorr, r3.Gaus(_track_Proj_true_t[iproj], sigmat), _track_Proj_true_t[iproj], isProbEScat );
            _track_TOFmeas[itrk].push_back(tempTOF);
          }
        } 
      }
      if (verbosityBASE > 3) std::cout << "timing layer count: " << _track_nTTL[itrk] << std::endl;
      if ( (_track_Proj_Clas[iproj] < 4 && _track_Proj_true_t[iproj] > 1e-3) ||  _track_Proj_Clas[iproj] == 5 || _track_Proj_Clas[iproj] == 6  )
        _track_RefProjID[itrk].push_back(iproj);
      
      if (_track_Proj_Clas[iproj] == 5 && (int)_track_trueID[itrk] > 0){
        projStrct tempProj;
        if ( IsECalProjection(_track_ProjLayer[iproj], _useAlternateForProjections) && trackSource == 0  && _track_Proj_t[iproj] >= 0  ){
          TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
          tempProj.caloid   = ReturnCalorimeterFromProjectionIndex(_track_ProjLayer[iproj]);
          tempProj.eta_Calo = projvec.Eta();
          tempProj.phi_Calo = projvec.Phi();
          _mcpart_EcalProjs[(int)_track_trueID[itrk]].push_back(tempProj);
        }

        if ( IsHCalProjection(_track_ProjLayer[iproj], _useAlternateForProjections) && trackSource == 0 && _track_Proj_t[iproj] >= 0 ){
          TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
          tempProj.caloid   = ReturnCalorimeterFromProjectionIndex(_track_ProjLayer[iproj]);
          tempProj.eta_Calo = projvec.Eta();
          tempProj.phi_Calo = projvec.Phi();
          _mcpart_HcalProjs[(int)_track_trueID[itrk]].push_back(tempProj);
        }
      }
      nCurrProj = iproj;      
    }
    if (trackSource == 1 )
      _track_hasOL[itrk] = (_track_hasOL[itrk] || true);
    if (trackSource == 2 )
      _track_hasIL[itrk]   = (_track_hasIL[itrk] || true);
    if ((int)_track_trueID[itrk] < 0) continue;
    if (_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size() > 1){
      if (verbosityBASE > 2) std::cout << "identified multiple tracks for: " << (int)_track_trueID[itrk] << " current track id: " << itrk << std::endl;
      for (int rtr = 0; rtr < (int)_mcpart_RecTrackIDs[(int)_track_trueID[itrk]].size(); rtr++){
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
    
    _track_matchECal[itrk] = matchingCalStrct();
    _track_matchHCal[itrk] = matchingCalStrct();
  }
  if (verbosityBASE > 0)
    if (trackID_ScatteredElectron != -1) std::cout << "found scattered electron: " << trackID_ScatteredElectron << std::endl;
  
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
    if (_mcpart_EcalProjs[imc].size() > 0){
       _mcpart_EcalProjs[imc].clear();
      _mcpart_EcalProjs[imc].resize(0);
    }
    if (_mcpart_HcalProjs[imc].size() > 0){
       _mcpart_HcalProjs[imc].clear();
      _mcpart_HcalProjs[imc].resize(0);
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
    if (_track_TOFmeas[itrk].size() > 0){
       _track_TOFmeas[itrk].clear();
      _track_TOFmeas[itrk].resize(0);
    }
    
  }
  for (Int_t iproj=0; iproj<_nProjections; iproj++){
    _track_Proj_Clas[iproj] = 0;
  }
}

