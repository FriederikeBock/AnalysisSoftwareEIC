#include <TLorentzVector.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

int verbosity_DISREC = 0;
int missedScatteredElectronDISLeptonReco = 0;
// Born information
double born_xB;
double born_Q2;
double born_y;
double born_W2;
double born_eta;

// HEPMC info
double hepmc_xB;
double hepmc_Q2;

// Born information
double lepton_xB;
double lepton_Q2;
double lepton_y;
double lepton_W2;
double lepton_eta;

// Reconstructed track information
int fRecElectron;
double e_tr_px;	
double e_tr_py;	
double e_tr_pz;	
double e_tr_p;	
double e_tr_th;
double e_tr_ph;
int e_tr_PDG;

// Reconstructed track information (require primary)
int fRecElectronPrimary;
double e_tr_px_primary;	
double e_tr_py_primary;	
double e_tr_pz_primary;	
double e_tr_p_primary;	
double e_tr_th_primary;
double e_tr_ph_primary;
            
// Electron reco
double e_tr_y;
double e_tr_Q2;
double e_tr_xB;
double e_tr_W2;
double e_tr_eta;

// Electron reco (require primary)
double e_tr_y_primary;
double e_tr_Q2_primary;
double e_tr_xB_primary;
double e_tr_W2_primary;
double e_tr_eta_primary;

// Hadronic quantities
double E_h;
double pz_h;
double delta_h;
double pT_h;
double gamma_h;
double transverse_E;
double missing_pT;
double E_tot;
double pz_tot;

// JB reco
double JB_y;
double JB_Q2;	
double JB_xB;	
double JB_W2;	
double JB_eta;

// JB2 reco
double JB2_y;
double JB2_Q2;	
double JB2_xB;	
double JB2_W2;	
double JB2_eta;

// DA reco
double DA_y;
double DA_Q2;	
double DA_xB;	
double DA_W2;	
double DA_eta;

// JB reco with tracks
double JB_tr_y;
double JB_tr_Q2;	
double JB_tr_xB;	
double JB_tr_W2;	
double JB_tr_eta;

// JB2 reco with tracks
double JB2_tr_y;
double JB2_tr_Q2;	
double JB2_tr_xB;	
double JB2_tr_W2;	
double JB2_tr_eta;

// DA reco with tracks
double DA_tr_y;
double DA_tr_Q2;	
double DA_tr_xB;	
double DA_tr_W2;	
double DA_tr_eta;

/*
// GEANT4 track information
double g4_tr_px;	
double g4_tr_py;	
double g4_tr_pz;	
double g4_tr_p;	
double g4_tr_th;	
double g4_tr_ph;	
	 
// GEANT4 DIS kinematics
double g4_tr_y;
double g4_tr_Q2;
double g4_tr_xB;
double g4_tr_W2;
double g4_tr_eta;
*/

double fEOP = 0.775;

const int NEMCALS = 3;
int EMCALIDS[NEMCALS] = {1, 3, 10};

double Me = 0.000510999;	// GeV
double Mp = 0.938272088;  	// GeV
double Mpi = 0.13957039;  	// GeV
double Mk = 0.493677;		// GeV

double fElectronE = 10.;
double fHadronE = 100.;

double fCrossingAngle = 25.e-3; 

int fTrackSource = 0;
int fClusterAlgo = 0;

TLorentzVector e_beam;
TLorentzVector h_beam;

TString fDISFilename;
TFile* disFile;
TTree* disTree;

TFile* caloResFile;
TH2D* caloRes;

TFile* trackResFile;
TH2D* trackRes;

bool useTruePDISrecPID = true;
bool use1DEoPDISrecPID = false;

TFile* rejFile; 
TDirectoryFile* eemcFile;
TDirectoryFile* becalFile;
TDirectoryFile* femcFile; 

TGraph* EopEEMC;
TGraph* EopBECAL;
TGraph* EopFEMC;
TH2D* eopcut2Dhist;

void SetBeamEnergies(double, double);
void SetCrossingAngle(double);
void InitializeDISRecon();
void InitializeDISTree();

void LoadEopCutGraph();
void LoadResolutions();

void ResetDISVariables();
bool LeptonReconstruction();
void HadronReconstruction();
void GetBornEvent();
int ElectronPID(bool);
void GetMissingTransverse();

void WriteDISTree( int nEventsTotal);

TLorentzVector CorrectCrossingAngle(TLorentzVector);
void CalculateElectronKinematics(TLorentzVector lp, double& xB, double& Q2, double& W2, double& y, double& eta);
void CalculateJacquetBlondelKinematics(TLorentzVector finalHadron, double& xB, double& Q2, double& W2, double& y, double& eta);
void CalculateJacquetBlondel2Kinematics(TLorentzVector finalHadron, double& xB, double& Q2, double& W2, double& y, double& eta);
void CalculateDoubleAngleKinematics(TLorentzVector lp, TLorentzVector finalHadron, double& xB, double& Q2, double& W2, double& y, double& eta);

