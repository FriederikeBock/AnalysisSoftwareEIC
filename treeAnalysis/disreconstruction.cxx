#include "disreconstruction.h"
#include <vector>

using namespace std;

void disreconstruction() {

  ResetDISVariables();

  GetBornEvent();
  bool havelepton =  LeptonReconstruction();
//	LeptonReconstructionPrimary();
//	HadronReconstruction();
  GetMissingTransverse();

  disTree->Fill();
}

void SetBeamEnergies(double electronE, double hadronE) {

  fElectronE = electronE;
  fHadronE = hadronE;

}

void SetCrossingAngle(double xingAngle) {
  fCrossingAngle = xingAngle;
}


void SetDISFilename(TString disfilename) {
  fDISFilename = disfilename;
}

void InitializeDISRecon() {

  std::cout << "Setting beam configuration to " << fElectronE << "x" << fHadronE << " GeV" << std::endl;
  e_beam = TLorentzVector(0.,0.,-sqrt(fElectronE * fElectronE - Me*Me), fElectronE);
  h_beam = TLorentzVector(0.,0.,sqrt(fHadronE*fHadronE - Mp*Mp), fHadronE);

  disFile = new TFile(fDISFilename, "RECREATE");
  InitializeDISTree();

  LoadEopCutGraph();
  LoadResolutions();
  
}

void LoadEopCutGraph() {

  // 1D
  rejFile   = new TFile("inputs/output_EoverP.root");
  eemcFile  = (TDirectoryFile*)rejFile->Get("EEMC");
  becalFile = (TDirectoryFile*)rejFile->Get("BECAL");
  femcFile  = (TDirectoryFile*)rejFile->Get("FEMC");

  TString eopCutType  = "95effi";
  TString pType       = "trueP";
  if (!useTruePDISrecPID)
    pType       = "recP";
  
  EopEEMC   = (TGraph*)eemcFile->Get(Form("graphEPcutValue_%s_%s_EEMC", pType.Data(), eopCutType.Data()));
  EopBECAL  = (TGraph*)becalFile->Get(Form("graphEPcutValue_%s_%s_BECAL", pType.Data(), eopCutType.Data()));
  EopFEMC   = (TGraph*)femcFile->Get(Form("graphEPcutValue_%s_%s_FEMC", pType.Data(), eopCutType.Data()));

  // 2D
  std::cout << Form("eopCut2D_%s_%s", pType.Data(), eopCutType.Data()) << std::endl;
  eopcut2Dhist  = (TH2D*)rejFile->Get(Form("eopCut2D_%s_%s", pType.Data(), eopCutType.Data()));
  if (!eopcut2Dhist) std::cout << "failed to find 2D hist" << std::endl;
}

void LoadResolutions() {

  caloResFile = new TFile("inputs/calo_resolution.root");
  caloRes = (TH2D*)caloResFile->Get("electron_resoFWHM");

  trackResFile = new TFile("inputs/track_resolution.root");
  trackRes = (TH2D*)trackResFile->Get("Electron_reso_track");
}

bool LeptonReconstruction() {

  int eTrackID = ElectronPID(false);	

  if( eTrackID < 0 ){
    if (verbosity_DISREC > 0) std::cout << "aborting Lepton reco" << std::endl;
    return false;
  }

  fRecElectron = 1;

  e_tr_PDG = (int)_mcpart_PDG[(int)_track_trueID[eTrackID]];
  
  TVector3 e_tr_vec(_track_px[eTrackID], _track_py[eTrackID], _track_pz[eTrackID]);
  double track_mag = e_tr_vec.Mag();
  double track_eta = e_tr_vec.PseudoRapidity();
  
  int caloID      = _track_matchECal[eTrackID].caloid;
  int clustID     = _track_matchECal[eTrackID].id;
  double calo_mag = (_clusters_calo[fClusterAlgo][caloID].at(clustID)).cluster_E;
  
  int track_res_bin = trackRes->FindBin(track_eta, track_mag); 
  int calo_res_bin = caloRes->FindBin(track_eta, track_mag); 
  
  double sigma_track = trackRes->GetBinContent(track_res_bin);
  double sigma_calo  =  caloRes->GetBinContent(calo_res_bin);
  
  double track_weight = 1. / (sigma_track * sigma_track);
  double calo_weight = 1. / (sigma_calo * sigma_calo);
  
  double average_mag = (track_mag*track_weight + calo_mag*calo_weight) / (track_weight + calo_weight);
  e_tr_vec.SetMag(average_mag);

  e_tr_p  = e_tr_vec.Mag();
  e_tr_th = e_tr_vec.Theta();
  e_tr_ph = e_tr_vec.Phi();

  e_tr_px = e_tr_vec.X();
  e_tr_py = e_tr_vec.Y();
  e_tr_pz = e_tr_vec.Z();
    
  TLorentzVector rec_eprime = CorrectCrossingAngle(TLorentzVector(e_tr_vec.x(), e_tr_vec.y(), e_tr_vec.z(), sqrt(e_tr_vec.Mag2() + Me*Me)));

  CalculateElectronKinematics(rec_eprime, e_tr_xB, e_tr_Q2, e_tr_W2, e_tr_y, e_tr_eta);	
  
  return true;

}

void LeptonReconstructionPrimary() {

  int eTrackID = ElectronPID(true);	
  if(eTrackID == -1) return;

  fRecElectronPrimary = 1;

  TVector3 e_tr_vec(_track_px[eTrackID], _track_py[eTrackID], _track_pz[eTrackID]);

  e_tr_p_primary  = e_tr_vec.Mag();
  e_tr_th_primary = e_tr_vec.Theta();
  e_tr_ph_primary = e_tr_vec.Phi();

  e_tr_px_primary = e_tr_vec.X();
  e_tr_py_primary = e_tr_vec.Y();
  e_tr_pz_primary = e_tr_vec.Z();
    
  TLorentzVector rec_eprime = CorrectCrossingAngle(TLorentzVector(e_tr_vec.x(), e_tr_vec.y(), e_tr_vec.z(), sqrt(e_tr_vec.Mag2() + Me*Me)));

  CalculateElectronKinematics(rec_eprime, e_tr_xB_primary, e_tr_Q2_primary, e_tr_W2_primary, e_tr_y_primary, e_tr_eta_primary);	

  return;

}


void HadronReconstruction() {

  int eTrackID = ElectronPID(false);

  TVector3 e_tr_vec(_track_px[eTrackID], _track_py[eTrackID], _track_pz[eTrackID]);
  TLorentzVector rec_eprime = CorrectCrossingAngle(TLorentzVector(e_tr_vec.x(), e_tr_vec.y(), e_tr_vec.z(), sqrt(e_tr_vec.Mag2() + Me*Me)));

  // pure JB/DA
  vector<TLorentzVector> clusters_all;

  // JB/DA with tracks when available
  vector<TLorentzVector> clusters_unmatched;
  vector<TLorentzVector> tracks;

  for(int icalo = 0; icalo < _active_calo; icalo++) {
    
    if(!caloEnabled[icalo]) continue;	

    for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[fClusterAlgo][icalo].size(); iclus++) {

      bool electronCluster = false;
      bool hasTrack = false;

      // Default pion mass
      double clusterMass = Mpi; 
      // But let's see if we can get PID from a matched track
      if((_clusters_calo[0][icalo].at(iclus)).cluster_isMatched) {

        hasTrack = true;
      
        for(int itrm = 0; itrm < (int)((_clusters_calo[0][icalo].at(iclus)).cluster_matchedTracks.size()); itrm++) {
          int matchedTrackID = (_clusters_calo[0][icalo].at(iclus)).cluster_matchedTracks[itrm].id;

          if(matchedTrackID == eTrackID) {
            electronCluster = true;
          }

          for(int itrk = 0; itrk < _nTracks; itrk++) {
              
            if(_track_source[itrk] != fTrackSource) continue;
              
            if(_track_ID[itrk] == matchedTrackID) {
              
              double pionLL = _track_pion_LL[itrk];
              double kaonLL = _track_kaon_LL[itrk];
              double protonLL = _track_proton_LL[itrk];

              if(pionLL > kaonLL && pionLL > protonLL && pionLL != -100 && pionLL < 0.) {
                clusterMass = Mpi;  	// it's a pion!
              } else if(kaonLL > pionLL && kaonLL > protonLL && kaonLL != 100 && kaonLL < 0.) {
                clusterMass = Mk;  	// it's a kaon!
              } else if(protonLL > pionLL && protonLL > kaonLL && protonLL != -100 && protonLL < 0.) {
                clusterMass = Mp;	// it's a proton!
              } 

            } // found our matching track
          } // track loop 
                                } // match track ID loop  
      } // is cluster matched						

      double clusterEta = (_clusters_calo[0][icalo].at(iclus)).cluster_Eta;
      double clusterTheta = 2. * atan(exp(-clusterEta));
      double clusterPhi = (_clusters_calo[0][icalo].at(iclus)).cluster_Phi;
          
      double clusterP = (_clusters_calo[0][icalo].at(iclus)).cluster_E;
      double clusterE = sqrt(clusterP*clusterP + clusterMass*clusterMass);

      double clusterPx = clusterP * sin(clusterTheta) * cos(clusterPhi);
      double clusterPy = clusterP * sin(clusterTheta) * sin(clusterPhi);
      double clusterPz = clusterP * cos(clusterTheta);
      
      // For pure JB...want all clusters except electron cluster
      if(!electronCluster) {
        clusters_all.push_back(CorrectCrossingAngle(TLorentzVector(clusterPx, clusterPy, clusterPz, clusterE)));
      }
      // For "track" JB...want all clusters not matched to a track
      if(!hasTrack) {
        clusters_unmatched.push_back(CorrectCrossingAngle(TLorentzVector(clusterPx, clusterPy, clusterPz, clusterE)));	
      }
    } // cluster loop
  } // calo loop	

  for(int itrk = 0; itrk < _nTracks; itrk++) {
      
    if(_track_source[itrk] != fTrackSource) continue;
        
    if(_track_ID[itrk] == eTrackID) continue;

    double trackMass = Mpi; 

    double pionLL = _track_pion_LL[itrk];
    double kaonLL = _track_kaon_LL[itrk];
    double protonLL = _track_proton_LL[itrk];

    if(pionLL > kaonLL && pionLL > protonLL && pionLL != -100 && pionLL < 0.) {
      trackMass = Mpi;  	// it's a pion!
    } else if(kaonLL > pionLL && kaonLL > protonLL && kaonLL != -100 && kaonLL < 0.) {
      trackMass = Mk;  	// it's a kaon!
    } else if(protonLL > pionLL && protonLL > kaonLL && protonLL != -100 && protonLL < 0.) {
      trackMass = Mp;		// it's a proton!
    } 

    double trackPx = _track_px[itrk];		
    double trackPy = _track_py[itrk];		
    double trackPz = _track_pz[itrk];		
    double trackP = sqrt(trackPx*trackPx + trackPy*trackPy + trackPz*trackPz);
    double trackE = sqrt(trackP*trackP + trackMass*trackMass);
    
    tracks.push_back(TLorentzVector(trackPx, trackPy, trackPz, trackE));

  } // track loop 

  // Sum clusters for pure JB/DA
  TLorentzVector sumHadron = TLorentzVector(0.,0.,0.,0.);
  for(int cl = 0; cl < (int)clusters_all.size(); cl++) {
    sumHadron+=clusters_all[cl];
  }

  // Sum clusters and tracks for JB/DA with tracks
  TLorentzVector sumHadronTr = TLorentzVector(0.,0.,0.,0.);
  for(int tr = 0; tr < (int)tracks.size(); tr++) {
    sumHadronTr+=tracks[tr];
  }
  for(int cl = 0; cl < (int)clusters_unmatched.size(); cl++) {
    sumHadronTr+=clusters_unmatched[cl];
  }

  CalculateJacquetBlondelKinematics(sumHadron, JB_xB, JB_Q2, JB_W2, JB_y, JB_eta);
  CalculateJacquetBlondelKinematics(sumHadronTr, JB_tr_xB, JB_tr_Q2, JB_tr_W2, JB_tr_y, JB_tr_eta);
  
  CalculateJacquetBlondel2Kinematics(sumHadron, JB2_xB, JB2_Q2, JB2_W2, JB2_y, JB2_eta);
  CalculateJacquetBlondel2Kinematics(sumHadronTr, JB2_tr_xB, JB2_tr_Q2, JB2_tr_W2, JB2_tr_y, JB2_tr_eta);

  if(!(eTrackID==-1)) {
    CalculateDoubleAngleKinematics(rec_eprime, sumHadron, DA_xB, DA_Q2, DA_W2, DA_y, DA_eta);
    CalculateDoubleAngleKinematics(rec_eprime, sumHadronTr, DA_tr_xB, DA_tr_Q2, DA_tr_W2, DA_tr_y, DA_tr_eta);
  }
  
  return;

}

void GetMissingTransverse() {

  vector<TLorentzVector> clusters_all;

  if (verbosity_DISREC) std::cout <<   "missing transverse" << std::endl;

  for(int icalo = 0; icalo < _active_calo; icalo++) {
    
    if(!caloEnabled[icalo]) continue;	

    for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[fClusterAlgo][icalo].size(); iclus++) {

      // Default pion mass
      double clusterMass = Mpi; 
      // But let's see if we can get PID from a matched track
      if((_clusters_calo[0][icalo].at(iclus)).cluster_isMatched) {

        for(Int_t itrm = 0; itrm < (Int_t)((_clusters_calo[0][icalo].at(iclus)).cluster_matchedTracks.size()); itrm++) {
          int matchedTrackID = (_clusters_calo[0][icalo].at(iclus)).cluster_matchedTracks[itrm].id;

          for(int itrk = 0; itrk < _nTracks; itrk++) {
              
            if(_track_source[itrk] != fTrackSource) continue;
              
            if(_track_ID[itrk] == matchedTrackID) {
              
              double pionLL = _track_pion_LL[itrk];
              double kaonLL = _track_kaon_LL[itrk];
              double protonLL = _track_proton_LL[itrk];

              if(pionLL > kaonLL && pionLL > protonLL && pionLL != -100 && pionLL < 0.) {
                clusterMass = Mpi;  	// it's a pion!
              } else if(kaonLL > pionLL && kaonLL > protonLL && kaonLL != 100 && kaonLL < 0.) {
                clusterMass = Mk;  	// it's a kaon!
              } else if(protonLL > pionLL && protonLL > kaonLL && protonLL != -100 && protonLL < 0.) {
                clusterMass = Mp;	// it's a proton!
              } 

            } // found our matching track
          } // track loop 
                                } // match track ID loop  
      } // is cluster matched						

      double clusterEta = (_clusters_calo[0][icalo].at(iclus)).cluster_Eta;
      double clusterTheta = 2. * atan(exp(-clusterEta));
      double clusterPhi = (_clusters_calo[0][icalo].at(iclus)).cluster_Phi;
          
      double clusterP = (_clusters_calo[0][icalo].at(iclus)).cluster_E;
      double clusterE = sqrt(clusterP*clusterP + clusterMass*clusterMass);

      double clusterPx = clusterP * sin(clusterTheta) * cos(clusterPhi);
      double clusterPy = clusterP * sin(clusterTheta) * sin(clusterPhi);
      double clusterPz = clusterP * cos(clusterTheta);
      
      clusters_all.push_back(CorrectCrossingAngle(TLorentzVector(clusterPx, clusterPy, clusterPz, clusterE)));
    } // cluster loop

  } // calo loop	

  TLorentzVector sum_all = TLorentzVector(0.,0.,0.,0.);;

  transverse_E = 0;

  for(int cl = 0; cl < (int)clusters_all.size(); cl++) {
    sum_all+=clusters_all[cl];
    transverse_E+=clusters_all[cl].Et();
  }
  missing_pT = sum_all.Pt();
  E_tot = sum_all.E();
  pz_tot = sum_all.Pz();

  if (verbosity_DISREC) std::cout <<   "missing transverse done" << std::endl;
  
  return;
}



int ElectronPID(bool reqParentID ) {
  
  double EopCut = 0.01;   // default initialization
  int eTrack = -1;

  vector<int> electron_candidates;
  vector<double> e_eta;
  vector<double> e_pT;
  vector<double> e_Eop;

  for(Int_t itrk = 0; itrk < _nTracks; itrk++) {

    // use track source = 0
    if(_track_source[itrk] != fTrackSource) continue;
    // only look at negative tracks
    if((int)_track_charge[itrk] > 0) continue;

/*
    int mcPartBCID = (int)_mcpart_BCID[(int)_track_trueID[itrk]];
    int parentID = -1;
    for(int imcp = 0; imcp < _nHepmcp; imcp++) {
      if(mcPartBCID == (int)_hepmcp_BCID[imcp]) {
        parentID = _hepmcp_m2[imcp];
      }
    } 
    if(reqParentID && (parentID != 10001) ) continue;
*/
    Int_t mcID = _track_trueID[itrk];
    if (mcID < 0) continue;
    double mc_PDG     = _mcpart_PDG[mcID] ;

    // Get track info
    TVector3 track_vec(_track_px[itrk], _track_py[itrk], _track_pz[itrk]);
    TLorentzVector mc_vec(_mcpart_px[mcID], _mcpart_py[mcID], _mcpart_pz[mcID], _mcpart_E[mcID]);

    double track_p    = track_vec.Mag();
    double track_pT   = track_vec.Pt();
    double track_eta  = track_vec.PseudoRapidity();
    double mc_p       = mc_vec.P();
    double mc_pT      = mc_vec.Pt();
    double mc_eta     = mc_vec.Eta();

    if (verbosity_DISREC) 
      std::cout << "Track = " << itrk << " track ID = " << _track_ID[itrk] << "\t rec P =" <<  track_p << "\t charge =" << (int)_track_charge[itrk] << " track true ID = " << mcID << "\t pdg = "<< mc_PDG << "\t true P = " << mc_p << std::endl; 
    
    int caloID  = _track_matchECal[itrk].caloid;
    if (caloID == -1) 
      continue;
    
    if (!use1DEoPDISrecPID){
      if (useTruePDISrecPID){
        EopCut = eopcut2Dhist->Interpolate(mc_eta, mc_p);
      } else {
        EopCut = eopcut2Dhist->Interpolate(mc_eta, track_p);
      }
    } else {
      if(caloID == kEEMC) {
        // EEMC
        if (useTruePDISrecPID)
          EopCut = max(0.01, EopEEMC->Eval(mc_p));
        else 
          EopCut = max(0.01, EopEEMC->Eval(track_p));
      } else if(caloID == kBECAL) {
        // BECAL
        if (useTruePDISrecPID)
          EopCut = max(0.01, EopBECAL->Eval(mc_p));
        else
          EopCut = max(0.01, EopBECAL->Eval(track_p));
      } else if(caloID == kFEMC) {
        // FEMC
        if (useTruePDISrecPID)
          EopCut = max(0.01, EopFEMC->Eval(mc_p));
        else
          EopCut = max(0.01, EopFEMC->Eval(track_p));
      }
    }
    int clustID = _track_matchECal[itrk].id;
    double Ecal_track = (_clusters_calo[fClusterAlgo][caloID].at(clustID)).cluster_E;
    double EoP = Ecal_track/track_p;
    if (useTruePDISrecPID)
      EoP = Ecal_track/mc_p;
    
    if (verbosity_DISREC > 1){
      std::cout << "\t\t E/p cut: " <<  EopCut << "\t E cls: " << Ecal_track << "\t E/p: " << EoP << "\t eta rec " << track_eta << "\t  eta mc " << mc_eta ;
    }

    if (EoP > EopCut) {
      electron_candidates.push_back(itrk);
      e_eta.push_back(track_eta);
      e_pT.push_back(track_pT);
      e_Eop.push_back(EoP);
      if (verbosity_DISREC > 1){
        std::cout << "\t passed" << std::endl;
        if (mc_PDG != 11)  std::cout << "======================\t this was not an electron \t ========================================" << std::endl;
        if (trackID_ScatteredElectron != -1 && trackID_ScatteredElectron == itrk)
          std::cout << "======================\t SUCCESS YOU FOUND the scattered electron \t ========================================" << std::endl;
      }
    } else {
      if (verbosity_DISREC > 1) {
        std::cout << "\t failed" << std::endl;
        if (trackID_ScatteredElectron != -1 && trackID_ScatteredElectron == itrk)
          std::cout << "======================\t YOU JUST MISSED the scattered electron \t ========================================" << std::endl;
      }
    }
  } // track loop
  
  

  double max_Eop = 0.;

  bool foundScatE = false;
  for(int iec = 0; iec < (int)electron_candidates.size(); iec++) {
    if(e_Eop[iec] > max_Eop) {
      max_Eop = e_Eop[iec];
      eTrack = electron_candidates[iec];
      if (trackID_ScatteredElectron != -1 && trackID_ScatteredElectron == eTrack)
        foundScatE = true;
    }
  }
  if (!foundScatE)
    missedScatteredElectronDISLeptonReco++;
  return eTrack;
}


void GetBornEvent() {

  // First get truth event info 
  // (should correspond to leptontic truth)
  hepmc_xB = _hepmcp_x2;
  hepmc_Q2 = _hepmcp_Q2;	

  // Next loop over truth particles to find virtual photon
  // (should correspond to hadronic/Born truth)
  for (int imcp = 0; imcp < _nHepmcp; imcp++) {
    
    // exchange boson
    if (_hepmcp_status[imcp] == 3 && _hepmcp_PDG[imcp] == 23) { 
      double pxg = _hepmcp_px[imcp];
      double pyg = _hepmcp_py[imcp];
      double pzg = _hepmcp_pz[imcp];
      double Eg  = _hepmcp_E[imcp];
      TLorentzVector virtual_photon;
      virtual_photon.SetPxPyPzE(pxg, pyg, pzg, Eg);
      CalculateElectronKinematics(e_beam - virtual_photon, born_xB, born_Q2, born_W2, born_y, born_eta);				
    }
    
    if(_hepmcp_status[imcp] == 1 && _hepmcp_PDG[imcp] == 11 && _hepmcp_m1[imcp] == 7) {
      double pxe = _hepmcp_px[imcp];
      double pye = _hepmcp_py[imcp];
      double pze = _hepmcp_pz[imcp];
      double Ee  = _hepmcp_E[imcp];
      TLorentzVector scattered_electron;
      scattered_electron.SetPxPyPzE(pxe, pye, pze, Ee);
      CalculateElectronKinematics(scattered_electron, lepton_xB, lepton_Q2, lepton_W2, lepton_y, lepton_eta);				
    }
  } // HepMC particle loop
}

void InitializeDISTree() {

  disTree = new TTree("T", "Tree with DIS kinematic variables");

  //	disTree->Branch("event",		&mEventCounter,	"event/I");

  disTree->Branch("rec_electron",	&fRecElectron,	"rec_electron/I");
  disTree->Branch("rec_electron_primary",	&fRecElectronPrimary,	"rec_electron_primary/I");

  disTree->Branch("hepmc_xB",	&hepmc_xB,	"hepmc_xB/D");
  disTree->Branch("hepmc_Q2",	&hepmc_Q2,	"hepmc_Q2/D");

  disTree->Branch("lepton_y",	&lepton_y,	"lepton_y/D");
  disTree->Branch("lepton_Q2",	&lepton_Q2,	"lepton_Q2/D");
  disTree->Branch("lepton_xB",	&lepton_xB,	"lepton_xB/D");
  disTree->Branch("lepton_W2",	&lepton_W2,	"lepton_W2/D");
  disTree->Branch("lepton_eta",	&lepton_eta,	"lepton_eta/D");

  disTree->Branch("born_y",	&born_y,	"born_y/D");
  disTree->Branch("born_Q2",	&born_Q2,	"born_Q2/D");
  disTree->Branch("born_xB",	&born_xB,	"born_xB/D");
  disTree->Branch("born_W2",	&born_W2,	"born_W2/D");
  disTree->Branch("born_eta",	&born_eta,	"born_eta/D");

  disTree->Branch("e_tr_px",	&e_tr_px,	"e_tr_px/D");	
  disTree->Branch("e_tr_py",	&e_tr_py,	"e_tr_py/D");	
  disTree->Branch("e_tr_pz",	&e_tr_pz,	"e_tr_pz/D");	
  disTree->Branch("e_tr_p",	&e_tr_p,	"e_tr_p/D");	
  disTree->Branch("e_tr_th",	&e_tr_th,	"e_tr_th/D");	
  disTree->Branch("e_tr_ph",	&e_tr_ph,	"e_tr_ph/D");	
  disTree->Branch("e_tr_PDG",	&e_tr_PDG,	"e_tr_PDG/I");

  disTree->Branch("e_tr_y",	&e_tr_y,	"e_tr_y/D");
  disTree->Branch("e_tr_Q2",	&e_tr_Q2,	"e_tr_Q2/D");
  disTree->Branch("e_tr_xB",	&e_tr_xB,	"e_tr_xB/D");
  disTree->Branch("e_tr_W2",	&e_tr_W2,	"e_tr_W2/D");
  disTree->Branch("e_tr_eta",	&e_tr_eta,	"e_tr_eta/D");

  disTree->Branch("e_tr_px_primary",	&e_tr_px_primary,	"e_tr_px_primary/D");	
  disTree->Branch("e_tr_py_primary",	&e_tr_py_primary,	"e_tr_py_primary/D");	
  disTree->Branch("e_tr_pz_primary",	&e_tr_pz_primary,	"e_tr_pz_primary/D");	
  disTree->Branch("e_tr_p_primary",	&e_tr_p_primary,	"e_tr_p_primary/D");	
  disTree->Branch("e_tr_th_primary",	&e_tr_th_primary,	"e_tr_th_primary/D");	
  disTree->Branch("e_tr_ph_primary",	&e_tr_ph_primary,	"e_tr_ph_primary/D");	

  disTree->Branch("e_tr_y_primary",	&e_tr_y_primary,	"e_tr_y_primary/D");
  disTree->Branch("e_tr_Q2_primary",	&e_tr_Q2_primary,	"e_tr_Q2_primary/D");
  disTree->Branch("e_tr_xB_primary",	&e_tr_xB_primary,	"e_tr_xB_primary/D");
  disTree->Branch("e_tr_W2_primary",	&e_tr_W2_primary,	"e_tr_W2_primary/D");
  disTree->Branch("e_tr_eta_primary",	&e_tr_eta_primary,	"e_tr_eta_primary/D");


  disTree->Branch("E_h",		&E_h,		"E_h/D");
  disTree->Branch("pz_h",		&pz_h,		"pz_h/D");
  disTree->Branch("delta_h",	&delta_h,	"delta_h/D");
  disTree->Branch("pT_h",		&pT_h,		"pT_h/D");
  disTree->Branch("gamma_h",	&gamma_h,	"gamma_h/D");
  disTree->Branch("transverse_E",	&transverse_E,	"transverse_E/D");
  disTree->Branch("missing_pT",	&missing_pT,	"missing_pT/D");
  disTree->Branch("E_tot",	&E_tot,		"E_tot/D");
  disTree->Branch("pz_tot",	&pz_tot,	"pz_tot/D");

  disTree->Branch("JB_y",		&JB_y,		"JB_y/D");
  disTree->Branch("JB_Q2",	&JB_Q2,		"JB_Q2/D");
  disTree->Branch("JB_xB",	&JB_xB,		"JB_xB/D");
  disTree->Branch("JB_W2",	&JB_W2,		"JB_W2/D");
  disTree->Branch("JB_eta",	&JB_eta,	"JB_eta/D");

  disTree->Branch("JB2_y",	&JB2_y,		"JB2_y/D");
  disTree->Branch("JB2_Q2",	&JB2_Q2,	"JB2_Q2/D");
  disTree->Branch("JB2_xB",	&JB2_xB,	"JB2_xB/D");
  disTree->Branch("JB2_W2",	&JB2_W2,	"JB2_W2/D");
  disTree->Branch("JB2_eta",	&JB2_eta,	"JB2_eta/D");

  disTree->Branch("DA_y",		&DA_y,		"DA_y/D");
  disTree->Branch("DA_Q2",	&DA_Q2,		"DA_Q2/D");
  disTree->Branch("DA_xB",	&DA_xB,		"DA_xB/D");
  disTree->Branch("DA_W2",	&DA_W2,		"DA_W2/D");
  disTree->Branch("DA_eta",	&DA_eta,	"DA_eta/D");

  disTree->Branch("JB_tr_y",	&JB_tr_y,	"JB_tr_y/D");
  disTree->Branch("JB_tr_Q2",	&JB_tr_Q2,	"JB_tr_Q2/D");
  disTree->Branch("JB_tr_xB",	&JB_tr_xB,	"JB_tr_xB/D");
  disTree->Branch("JB_tr_W2",	&JB_tr_W2,	"JB_tr_W2/D");
  disTree->Branch("JB_tr_eta",	&JB_tr_eta,	"JB_tr_eta/D");

  disTree->Branch("JB2_tr_y",	&JB2_tr_y,	"JB2_tr_y/D");
  disTree->Branch("JB2_tr_Q2",	&JB2_tr_Q2,	"JB2_tr_Q2/D");
  disTree->Branch("JB2_tr_xB",	&JB2_tr_xB,	"JB2_tr_xB/D");
  disTree->Branch("JB2_tr_W2",	&JB2_tr_W2,	"JB2_tr_W2/D");
  disTree->Branch("JB2_tr_eta",	&JB2_tr_eta,	"JB2_tr_eta/D");

  disTree->Branch("DA_tr_y",	&DA_tr_y,	"DA_tr_y/D");
  disTree->Branch("DA_tr_Q2",	&DA_tr_Q2,	"DA_tr_Q2/D");
  disTree->Branch("DA_tr_xB",	&DA_tr_xB,	"DA_tr_xB/D");
  disTree->Branch("DA_tr_W2",	&DA_tr_W2,	"DA_tr_W2/D");
  disTree->Branch("DA_tr_eta",	&DA_tr_eta,	"DA_tr_eta/D");
}


void ResetDISVariables() {

  fRecElectron = 0;
  fRecElectronPrimary = 0;

  hepmc_xB = -99.;
  hepmc_Q2 = -99.;

  lepton_y = -99.; 
  lepton_Q2 = -99.; 
  lepton_xB = -99.; 
  lepton_W2 = -99.; 
  lepton_eta = -99.;

  born_y = -99.; 
  born_Q2 = -99.; 
  born_xB = -99.; 
  born_W2 = -99.; 
  born_eta = -99.;

  e_tr_px = -99.;
  e_tr_py = -99.;
  e_tr_pz = -99.;
  e_tr_p = -99.;
  e_tr_th = -99.;
  e_tr_ph = -99.;
  e_tr_PDG = -99.;
        
  e_tr_y = -99.;
  e_tr_Q2 = -99.;
  e_tr_xB = -99.;
  e_tr_W2 = -99.;
  e_tr_eta = -99.;

  e_tr_px_primary = -99.;
  e_tr_py_primary = -99.;
  e_tr_pz_primary = -99.;
  e_tr_p_primary = -99.;
  e_tr_th_primary = -99.;
  e_tr_ph_primary = -99.;
        
  e_tr_y_primary = -99.;
  e_tr_Q2_primary = -99.;
  e_tr_xB_primary = -99.;
  e_tr_W2_primary = -99.;
  e_tr_eta_primary = -99.;

  E_h = -99.;
  pz_h = -99.;
  delta_h = -99.;
  pT_h = -99.;
  gamma_h = -99.;
  transverse_E = -99.;
  missing_pT = -99.;
  E_tot = -99.;
  pz_tot = -99.;

  JB_y = -99.;
  JB_Q2 = -99.;  
  JB_xB = -99.;  
  JB_W2 = -99.;  
  JB_eta = -99.;

  JB2_y = -99.;
  JB2_Q2 = -99.;  
  JB2_xB = -99.;  
  JB2_W2 = -99.;  
  JB2_eta = -99.;

  DA_y = -99.;
  DA_Q2 = -99.;  
  DA_xB = -99.;  
  DA_W2 = -99.;  
  DA_eta = -99.;

  JB_tr_y = -99.;
  JB_tr_Q2 = -99.;  
  JB_tr_xB = -99.;  
  JB_tr_W2 = -99.;  
  JB_tr_eta = -99.;

  JB2_tr_y = -99.;
  JB2_tr_Q2 = -99.;  
  JB2_tr_xB = -99.;  
  JB2_tr_W2 = -99.;  
  JB2_tr_eta = -99.;

  DA_tr_y = -99.;
  DA_tr_Q2 = -99.;  
  DA_tr_xB = -99.;  
  DA_tr_W2 = -99.;  
  DA_tr_eta = -99.;


}

void WriteDISTree( int nEventsTotal ) {

  disFile->cd();
  disTree->Write("T", TObject::kOverwrite);
  disFile->Close();

  std::cout << "Fraction of missed scattered electron w/ calo PID: " << (float)(missedScatteredElectronDISLeptonReco)/(float)nEventsTotal << std::endl;
  
}

TLorentzVector CorrectCrossingAngle(TLorentzVector vecIn) {
  TLorentzVector vecOut = vecIn;
  vecOut.RotateY(fCrossingAngle/2.);
  vecOut.Boost(sin(fCrossingAngle/2.),0,0);
  return vecOut;
}

void CalculateElectronKinematics(TLorentzVector lp, double& xB, double& Q2, double& W2, double& y, double& eta) {

  TLorentzVector q = e_beam - lp;
  Q2 = -q*q;
  xB = Q2/(2. * (h_beam*q));
  W2 = (h_beam+q)*(h_beam+q);
  y = (q*h_beam)/(e_beam*h_beam);
  eta = lp.PseudoRapidity();

}


void CalculateJacquetBlondelKinematics(TLorentzVector finalHadron, double& xB, double& Q2, double& W2, double& y, double& eta) {

  // Now calculate variables
  y = ( finalHadron.E() - finalHadron.Pz() ) / ( 2.*e_beam.E() );
  Q2 = finalHadron.Perp2() / (1. - JB_y);
  double s = (e_beam + h_beam).Mag2();
  xB = Q2/(s*y);
  W2 = (Mp*Mp) + (Q2/xB) - Q2;
  double etheta = acos(1. - Q2/(2.*e_beam.E()*e_beam.E()*(1. - y)));

  // NOTE! No minus sign because
  // electron direction is negative eta
  eta = log(tan(etheta/2.));
    
}

void CalculateJacquetBlondel2Kinematics(TLorentzVector finalHadron, double& xB, double& Q2, double& W2, double& y, double& eta) {

  // Now calculate variables
  y = ( h_beam.E() * finalHadron.E() - h_beam.Pz() * finalHadron.Pz() - Mp*Mp) / ( e_beam.E() * h_beam.E() - e_beam.Pz() * h_beam.Pz() );
  Q2 = finalHadron.Perp2() / (1. - JB_y);
  double s = (e_beam + h_beam).Mag2();
  xB = Q2/(s*y);
  W2 = (Mp*Mp) + (Q2/xB) - Q2;
  double etheta = acos(1. - Q2/(2.*e_beam.E()*e_beam.E()*(1. - y)));

  // NOTE! No minus sign because
  // electron direction is negative eta
  eta = log(tan(etheta/2.));

}



void CalculateDoubleAngleKinematics(TLorentzVector lp, TLorentzVector finalHadron, double& xB, double& Q2, double& W2, double& y, double& eta) {

  E_h = finalHadron.E();
  pz_h = finalHadron.Pz();
  delta_h = E_h - pz_h;
  pT_h = finalHadron.Pt();
  gamma_h = 2. * atan(delta_h/pT_h);

  double theta = lp.Theta();
  
  y = tan(gamma_h/2.) / (tan(theta/2.) + tan(gamma_h/2.));
  Q2 = 4. * e_beam.E()*e_beam.E() / (tan(theta/2.) * (tan(theta/2.) + tan(gamma_h/2.)));
  double s = (e_beam + h_beam).Mag2();
  xB = Q2/(s*y);
  W2 = (Mp*Mp) + (Q2/xB) - Q2;	
  double etheta = acos(1. - Q2/(2.*e_beam.E()*e_beam.E()*(1. - y)));

  eta = lp.PseudoRapidity();

}


