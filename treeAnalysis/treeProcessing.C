#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include "treeProcessing.h"
#include "event_utils.cxx"
#include "jet_finder.cxx"
#include "jet_observables.cxx"
#include "caloheader.h"
#include "clusterizer.cxx"

#include "jetresolutionhistos.cxx"
#include "caloresolutionhistos.cxx"
#include "clusterstudies.cxx"
#include "trackingefficiency.cxx"
#include "hitstudies.cxx"
#include "trackmatchingstudies.cxx"

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <algorithm>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

void treeProcessing(
    TString inFile              = "",
    TString inFileGeometry      = "geometry.root",
    TString addOutputName       = "",
    std::string baseOutputDir   = ".",
    bool do_reclus              = true,
    bool do_jetfinding          = false,
    // Double_t maxNEvent = 1e5,
    bool hasTiming              = true,
    bool isALLSILICON           = true,
    Int_t maxNEvent             = -1,
    Int_t verbosity             = 0,
    bool doCalibration          = false,
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0,
    std::string jetAlgorithm    = "anti-kt",
    std::vector<double> jetRParameters = {0.3, 0.5, 0.8, 1.0},
    double tracked_jet_max_pT   = 30,
    bool runCaloRes             = true
){
    // make output directory
    TString dateForOutput = ReturnDateStr();
    outputDir = Form("%s/treeProcessing/%s", baseOutputDir.c_str(), addOutputName.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir + "/etaphi");

    if(do_jetfinding) _do_jetfinding = true;

    // load tree
    TChain *const tt_event = new TChain("event_tree");
    if (inFile.EndsWith(".root")) {                     // are we loading a single root tree?
        std::cout << "loading a single root file" << std::endl;
        tt_event->AddFile(inFile);
    }
    else {                                              // or are we loading a bunch?
        std::cout << "loading a list of files at " << inFile << std::endl;
        std::ifstream files(inFile);
        std::string filePath;

        while (std::getline(files, filePath)) {
            std::cout << filePath << "\n";
            tt_event->AddFile(filePath.c_str());
        }
        files.close();
    }
    if(!tt_event){ std::cout << "tree not found... returning!"<< std::endl; return;}

    // // load geometry tree
    tt_geometry =  (TTree *) (new TFile(inFileGeometry.Data(), "READ"))->Get("geometry_tree");
    if(!tt_geometry){ cout << "geometry tree not found... returning!"<< endl; return;}
    // load all branches (see header)
    SetBranchAddressesTree(tt_event);
    SetBranchAddressesGeometryTree(tt_geometry);
    SetGeometryIndices();

    for (Int_t c = 0; c < 12; c++){
      cout << str_calorimeter[c] << "\t" << caloEnabled[c] << endl;
    }

    Long64_t nEntriesTree                 = tt_event->GetEntries();
    std::cout << "Number of events in tree: " << nEntriesTree << std::endl;
    if(maxNEvent>0 && maxNEvent<nEntriesTree){
        nEntriesTree = maxNEvent;
        std::cout << "Will only analyze first " << maxNEvent << " events in the tree..." << std::endl;
    }
    _do_TimingStudies         = hasTiming;
    _is_ALLSILICON            = isALLSILICON;
    _doClusterECalibration    = doCalibration;
    if(_doClusterECalibration){
        std::cout << "clusters will be energy-corrected and subsequently smeared to meet testbeam constant term!" << std::endl;
    }
    // Additional setup
    // Jet energy scale + resolution only support one R, so we select it here
    double jetRForJES = 0.5;
    if (std::find(jetRParameters.begin(), jetRParameters.end(), jetRForJES) == jetRParameters.end()) {
      std::cout << "Requested R=" << jetRForJES << " for JES/JER, but it's not in the jetRParameters. Adjust one of the settings.\n";
    }
    // Base setup for jet related observables (defined here since we need to access them everywhere)
    auto eventObservables = EventObservables();
    std::map<std::string, JetObservables> jetObservables;
    std::vector<std::string> nPDFNames = {"ep", "EPPS16nlo_CT14nlo_Au197"};
    TRandom3 eASelector(0);
    if (_do_jetfinding) {
      // Create jet level observables
      std::cout << "About to create jet observables.\n" << std::endl;
      for (auto pdfLabel : nPDFNames) {
        std::string fullLabel = (pdfLabel != "") ? ("_" + pdfLabel) : "";
        jetObservables.emplace("true" + fullLabel, JetObservables(JetType_t::full, pdfLabel, "true"));
        jetObservables.emplace("true_charged" + fullLabel, JetObservables(JetType_t::charged, pdfLabel, "true"));
        jetObservables.emplace("charged" + fullLabel, JetObservables(JetType_t::charged, pdfLabel));
        jetObservables.emplace("calo" + fullLabel, JetObservables(JetType_t::calo, pdfLabel));
        jetObservables.emplace("full" + fullLabel, JetObservables(JetType_t::full, pdfLabel));
      }
      // Initialization
      // Event level
      eventObservables.Init();

      // Jet level
      for (auto && [k, v] : jetObservables) {
        v.Init(jetRParameters);
      }
    }

    _nEventsTree=0;
    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {
        // load current event
        tt_event->GetEntry(i);
        _nEventsTree++;

        // processing progress info
        if(i>0 && nEntriesTree>100 && i%(nEntriesTree/(50))==0) std::cout << "//processed " << 100*(i)/nEntriesTree << "%"  << std::endl;
        if(verbosity>0){
          std::cout << "***********************************************************************************************************" << std::endl;
          std::cout << "event " << i << std::endl;
          std::cout << "***********************************************************************************************************" << std::endl;
        }
        // calculate useful quantities
        for(Int_t imc=0; imc<_nMCPart; imc++){
            TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
            _mcpart_Eta[imc]=truevec.Eta();
            _mcpart_Phi[imc]=truevec.Phi();
            int motherid = 0;
            for(Int_t imc2=0; imc2<_nMCPart; imc2++){
                if(_mcpart_ID[imc2]==_mcpart_ID_parent[imc]){
                    motherid=imc2;
                }
            }
            // if(_mcpart_PDG[imc]==111) cout << "pi0 found (" << _mcpart_ID[imc] << ") with mother: " << _mcpart_PDG[motherid] << " (" <<  _mcpart_ID_parent[imc] << ")" << endl;
        }
        if (verbosity > 0){
          for(int imc=0; imc<_nMCPart; imc++){
            std::cout << "MC part: \t" << imc << "\t E = " << _mcpart_E[imc] << "\tEta: " << _mcpart_Eta[imc]<< "\tPhi: " << _mcpart_Phi[imc]<<  std::endl;
          }
          for (Int_t icalo = 0; icalo < maxcalo; icalo++){
            if (caloEnabled[icalo])
              std::cout << "towers " <<   str_calorimeter[icalo].Data() << "\t" << ReturnMaxTowerCalo(icalo) << std::endl;
          }
        }

        Int_t nTowers[maxcalo] = {_nTowers_FHCAL, _nTowers_FEMC, _nTowers_DRCALO, _nTowers_EEMC, _nTowers_CEMC,
                                  _nTowers_EHCAL, _nTowers_HCALIN, _nTowers_HCALOUT, _nTowers_LFHCAL, _nTowers_EEMCG,
                                  _nTowers_BECAL  };

        // run clusterizers normal calos
        for (int cal = 0; cal < maxcalo; cal++){
          if(do_reclus && nTowers[cal] > 0 && caloEnabled[cal]){
            for (int algo = 0; algo < _active_algo; algo++){
              if(verbosity>1) cout << "clusterizing " << str_clusterizer[algo].Data() <<  " for " << str_calorimeter[cal].Data() << endl;
              runclusterizer(algo, cal, seedE[cal], aggE[cal], primaryTrackSource);
            }
          }
        }

        if(do_reclus && _nTowers_FOCAL && caloEnabled[kFOCAL]){
            float seed_E_FOCAL = 0.1;
            float aggregation_E_FOCAL = 0.001;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for FOCAL" << endl;
                runclusterizer(kMA, kFOCAL,seed_E_FOCAL, aggregation_E_FOCAL, primaryTrackSource);
            }
        }

        if(do_reclus && _nTowers_DRCALO && caloEnabled[kDRCALO]){ //do_V1clusterizerDRCALO
          if(verbosity>1) cout << "clusterizing V1 for DRCALO" << endl;
          runclusterizer(kV1, kDRCALO, seedE[kDRCALO], aggE[kDRCALO], primaryTrackSource);
        }
        if((do_reclus) && verbosity>1) cout << "done with clusterization!" << endl;

        // ANCHOR Hits loop variables:
        // float* _hits_x[ihit]
        // float* _hits_y[ihit]
        // float* _hits_z[ihit]
        // float* _hits_t[ihit]
        // float* _hits_edep[ihit]
        // for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
        //     if(verbosity>1) std::cout << "\tHIT: hit " << ihit << "\tin layer " << _hits_layerID[ihit] << "\twith X = " << _hits_x[ihit] << " cm" << std::endl;
        // }
        if(tracksEnabled) hitstudies(primaryTrackSource);

        // ANCHOR Track loop variables:
        // float* _track_ID[itrk]
        // float* _track_trueID[itrk]
        // float* _track_px[itrk]
        // float* _track_py[itrk]
        // float* _track_pz[itrk]
        std::vector<float> jetf_track_px;
        std::vector<float> jetf_track_py;
        std::vector<float> jetf_track_pz;
        std::vector<float> jetf_track_E;
        std::vector<float> jetf_full_px;
        std::vector<float> jetf_full_py;
        std::vector<float> jetf_full_pz;
        std::vector<float> jetf_full_E;
        std::vector<float> jetf_all_px;
        std::vector<float> jetf_all_py;
        std::vector<float> jetf_all_pz;
        std::vector<float> jetf_all_E;
        double massHypothesis = 0.139;
        //std::cout << "Tracks\n";
        for(Int_t itrk=0; itrk<_nTracks; itrk++){
            _track_trueID[itrk] = GetCorrectMCArrayEntry(_track_trueID[itrk]);
            // if(verbosity>1) std::cout << "\tTrack: track " << itrk << "\twith true ID " << _track_trueID[itrk] << "\tand X = " << _track_px[itrk] << " cm" << std::endl;

            // Skip tracks that aren't selected.
            if(_track_source[itrk] != primaryTrackSource) {
                continue;
            }

            if(_do_jetfinding){
                TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
                float Etrack = TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)+TMath::Power(_track_pz[itrk],2) + TMath::Power(massHypothesis, 2));
                //std::cout << itrk << "-> px=" << _track_px[itrk]
                //          << ", py=" << _track_py[itrk]
                //          << ", pz=" << _track_pz[itrk]
                //          << ", E=" << Etrack
                //          << ", trueID=" << _track_trueID[itrk]
                //          << "\n";
                // create track vector for jet finder
                float track_pT = Etrack / cosh(trackvec.Eta());
                if (track_pT < tracked_jet_max_pT) {
                    jetf_track_px.push_back(_track_px[itrk]);
                    jetf_track_py.push_back(_track_py[itrk]);
                    jetf_track_pz.push_back(_track_pz[itrk]);
                    jetf_track_E.push_back(Etrack);
                }
                jetf_full_px.push_back(_track_px[itrk]);
                jetf_full_py.push_back(_track_py[itrk]);
                jetf_full_pz.push_back(_track_pz[itrk]);
                jetf_full_E.push_back(Etrack);
                jetf_all_px.push_back(_track_px[itrk]);
                jetf_all_py.push_back(_track_py[itrk]);
                jetf_all_pz.push_back(_track_pz[itrk]);
                jetf_all_E.push_back(Etrack);
            }
        }
        // obtain labels for different track sources
        prepareMCMatchInfo();


        // ANCHOR jet variables:
        std::vector<float> jetf_hcal_E;
        std::vector<float> jetf_hcal_px;
        std::vector<float> jetf_hcal_py;
        std::vector<float> jetf_hcal_pz;
        std::vector<float> jetf_emcal_E;
        std::vector<float> jetf_emcal_px;
        std::vector<float> jetf_emcal_py;
        std::vector<float> jetf_emcal_pz;
        std::vector<float> jetf_calo_E;
        std::vector<float> jetf_calo_px;
        std::vector<float> jetf_calo_py;
        std::vector<float> jetf_calo_pz;

        int fwdCaloID = kFHCAL;
        if (!caloEnabled[kFHCAL]) fwdCaloID = kLFHCAL;

        for (Int_t icalo = 0; icalo < maxcalo; icalo++){
          for (Int_t ialgo = 0; ialgo < maxAlgo; ialgo++){
            if (loadClusterizerInput(ialgo, icalo) && verbosity>2) std::cout << str_calorimeter[icalo].Data() << "\t" << ialgo << endl;
          }
        }

        if (verbosity > 4){
          for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[kV1][fwdCaloID].size(); iclus++){
            std::cout << "\tcls " << iclus << "\tE " << (_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_E << "\tEta " << (_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_Eta << "\tPhi " << (_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_Phi << "\tntowers: " << (_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_NTowers << "\ttrueID: " << _clusters_FHCAL_trueID[iclus] << std::endl;
          }
        }

        // apply calibration if desired
        if(kV1<_active_algo){
          // TODO: Should kV1 in the calibration actually be _active_algo?
          for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[kV1][fwdCaloID].size(); iclus++){
            (_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_trueID = GetCorrectMCArrayEntry((_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_trueID);
            if(_doClusterECalibration){
              (_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_E/=getCalibrationValue((_clusters_calo[kV1][fwdCaloID].at(iclus)).cluster_E, fwdCaloID, kV1);
            }
          }
          for(Int_t iclus=0; iclus<(Int_t)_clusters_calo[kV1][kFEMC].size(); iclus++){
            (_clusters_calo[kV1][kFEMC].at(iclus)).cluster_trueID = GetCorrectMCArrayEntry((_clusters_calo[kV1][kFEMC].at(iclus)).cluster_trueID);
            if(_doClusterECalibration){
              (_clusters_calo[kV1][kFEMC].at(iclus)).cluster_E/=getCalibrationValue((_clusters_calo[kV1][kFEMC].at(iclus)).cluster_E, kFEMC, kV1);
            }
          }
        }
        if(do_reclus && kMA<_active_algo && _do_jetfinding){
          // ANCHOR FEMC cluster loop variables:
          // for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
          //     if(verbosity>1) std::cout << "\tFEMC:  cluster " << iclus << "\twith E = " << (_clusters_calo[kV1][kFEMC].at(iclus)).cluster_E << " GeV" << std::endl;
          // }
          if (caloEnabled[kFEMC]){
            fillECalClustersIntoJetFindingInputs(
                kFEMC, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
          }
          // EEMC
          if (caloEnabled[kEEMC]){
            fillECalClustersIntoJetFindingInputs(
                kEEMC, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
          }
          // EEMCG
          if (caloEnabled[kEEMCG]){
            fillECalClustersIntoJetFindingInputs(
                kEEMCG, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
          }
          // CEMC
          if (caloEnabled[kCEMC]){
            fillECalClustersIntoJetFindingInputs(
                kCEMC, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
          }
          // BECAL
          if (caloEnabled[kBECAL]){
            fillECalClustersIntoJetFindingInputs(
                kBECAL, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
          }

          // Forward HCal
          if (caloEnabled[fwdCaloID]){
            fillHCalClustersIntoJetFindingInputs(
                fwdCaloID, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
          }
          // EHCAL
          if (caloEnabled[kEHCAL]){
            fillHCalClustersIntoJetFindingInputs(
                kEHCAL, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
          }
          // HCALIN
          if (caloEnabled[kHCALIN]){
            fillHCalClustersIntoJetFindingInputs(
                kHCALIN, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
          }
          // HCALOUT
          if (caloEnabled[kHCALOUT]){
            fillHCalClustersIntoJetFindingInputs(
                kHCALOUT, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
          }
        }

        // ANCHOR MC particle loop variables:
        // float * _mcpart_ID[imc]
        // float * _mcpart_ID_parent[imc]
        // float * _mcpart_PDG[imc]
        // float * _mcpart_E [imc]
        // float * _mcpart_px[imc]
        // float * _mcpart_py[imc]
        // float * _mcpart_pz[imc]
        if(_do_jetfinding){
            std::vector<float> jetf_truth_px;
            std::vector<float> jetf_truth_py;
            std::vector<float> jetf_truth_pz;
            std::vector<float> jetf_truth_E;
            std::vector<float> jetf_truthcharged_px;
            std::vector<float> jetf_truthcharged_py;
            std::vector<float> jetf_truthcharged_pz;
            std::vector<float> jetf_truthcharged_E;
            for(Int_t imc=0; imc<_nMCPart; imc++){
                //std::cout << imc << "-> id=" << _mcpart_PDG[imc] << ", barcode=" << _mcpart_BCID[imc] << "\n";
                TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
//                 if(truevec.Eta()<1){
                    if(verbosity>1) std::cout << "\tMC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << std::endl;
                    //  create truth vector for jet finder
                    jetf_truth_px.push_back(_mcpart_px[imc]);
                    jetf_truth_py.push_back(_mcpart_py[imc]);
                    jetf_truth_pz.push_back(_mcpart_pz[imc]);
                    jetf_truth_E.push_back(_mcpart_E[imc]);
//                 }
                // fill charged only inputs (reject neutral particles)
//                 std::cout << _mcpart_PDG[imc] << std::endl;
                if(abs(_mcpart_PDG[imc])==211 || abs(_mcpart_PDG[imc])==321 || abs(_mcpart_PDG[imc])==2212 || abs(_mcpart_PDG[imc])==11 || abs(_mcpart_PDG[imc])==13){
                    jetf_truthcharged_px.push_back(_mcpart_px[imc]);
                    jetf_truthcharged_py.push_back(_mcpart_py[imc]);
                    jetf_truthcharged_pz.push_back(_mcpart_pz[imc]);
                    jetf_truthcharged_E.push_back(_mcpart_E[imc]);
                }
            }
            // Compare to hepmc
            //std::cout << "HepMC\n";
            //for (unsigned int i = 0; i < _nHepmcp; ++i) {
            //  std::cout << i << "-> id=" << _hepmcp_PDG[i] << ", barcode=" << _hepmcp_BCID[i] << ", px=" << _hepmcp_px[i]
            //            << ", py=" << _hepmcp_py[i]
            //            << ", pz=" << _hepmcp_pz[i]
            //            << ", E=" << _hepmcp_E[i]
            //            << ", status=" << _hepmcp_status[i]
            //            << "\n";
            //}
            //std::cout << "End\n";

            std::vector<double> nocluster_px;
            std::vector<double> nocluster_py;
            std::vector<double> nocluster_pz;
            std::vector<double> nocluster_E;
            if (caloEnabled[kFEMC]){
                for (uint32_t i = 0; i < (uint32_t)_nTowers_FEMC; i++) {
                    if (_tower_FEMC_E[i] < seedE[kFEMC]) {
                        continue;
                    }
                    TVector3 twrPos = TowerPositionVectorFromIndicesGeometry(_tower_FEMC_iEta[i], _tower_FEMC_iPhi[i], 0, kFEMC);
                    double pt = _tower_FEMC_E[i] / cosh(twrPos.Eta());
                    nocluster_px.push_back(pt * cos(twrPos.Phi()));
                    nocluster_py.push_back(pt * sin(twrPos.Phi()));
                    nocluster_pz.push_back(pt * sinh(twrPos.Eta()));
                    nocluster_E.push_back(_tower_FEMC_E[i]);
                }
            }
            // Event level
            auto kinematics = fillEventObservables(eventObservables, primaryTrackSource);
            auto hepmcCalculatedStored = KinematicsUsingTrueInfo(DirectKinematicsOptions_t::kHepMCStored, primaryTrackSource, verbosity);
            auto hepmcCalculatedKinematcs = KinematicsUsingTrueInfo(DirectKinematicsOptions_t::kHepMCCalculated, primaryTrackSource, verbosity);
            auto partLevelCalculatedKinematcs = KinematicsUsingTrueInfo(DirectKinematicsOptions_t::kPartLevel, primaryTrackSource, verbosity);
            try {
              auto detLevelCalculatedKinematcs = KinematicsUsingTrueInfo(DirectKinematicsOptions_t::kDetLevel, primaryTrackSource, verbosity);
            }
            catch (KinematicsErrors_t e) {
              std::cout << "Reco level kinematics error: " << e << ". Skipping...\n";
            }

            bool skipJetFindingDueToKinematics = false;
            if (kinematics.x < 0 || kinematics.x > 1) {
              skipJetFindingDueToKinematics = true;
            }

            // ANCHOR JET FINDING
            if (skipJetFindingDueToKinematics == false) {
              for (double jetR : jetRParameters) {
                // truth jets
                auto jetsTrue = findJets(jetR, jetAlgorithm, jetf_truth_px, jetf_truth_py, jetf_truth_pz, jetf_truth_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsTrue).size() << " true jets" << std::endl;        // printJets(std::get<1>(jetsTrue));

                auto jetsTrueCharged = findJets(jetR, jetAlgorithm, jetf_truthcharged_px, jetf_truthcharged_py, jetf_truthcharged_pz, jetf_truthcharged_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsTrueCharged).size() << " true charged jets" << std::endl;        // printJets(std::get<1>(jetsTrue));

                // track-based jets (rec)
                auto jetsTrackRec = findJets(jetR, jetAlgorithm, jetf_track_px, jetf_track_py, jetf_track_pz, jetf_track_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsTrackRec).size() << " rec track jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

                // full jets (rec)
                auto jetsFullRec = findJets(jetR, jetAlgorithm, jetf_full_px, jetf_full_py, jetf_full_pz, jetf_full_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsFullRec).size() << " rec full jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

                // hcal jets (rec)
                auto jetsHcalRec = findJets(jetR, jetAlgorithm, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz, jetf_hcal_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsHcalRec).size() << " rec hcal jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

                // emcal jets (rec)
                auto jetsEmcalRec = findJets(jetR, jetAlgorithm, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz, jetf_emcal_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsEmcalRec).size() << " rec emcal jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

                // calo jets (rec)
                auto jetsCaloRec = findJets(jetR, jetAlgorithm, jetf_calo_px, jetf_calo_py, jetf_calo_pz, jetf_calo_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsCaloRec).size() << " rec calo jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

                // all jets (rec)
                auto jetsAllRec = findJets(jetR, jetAlgorithm, jetf_all_px, jetf_all_py, jetf_all_pz, jetf_all_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsAllRec).size() << " rec all jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));

                // calo jets, no clustering
                auto jetsNoCluster = findJets(jetR, jetAlgorithm, nocluster_px, nocluster_py, nocluster_pz, nocluster_E);
                if(verbosity>1) std::cout << "found " << std::get<1>(jetsNoCluster).size() << " rec calo jets" << std::endl;        // printJets(std::get<1>(jetsTrackRec));


                // Jet observables
                // We want the ep and eA samples to be independent, so divide it here
                bool fill_eA = (eASelector.Rndm() >= 0.5);
                for (auto pdfLabel : nPDFNames) {
                  // IF we're filling eA, skip ep
                  if (fill_eA == true && pdfLabel == "ep") {
                    continue;
                  }
                  // If we're filling ep, skip eA
                  if (fill_eA == false && pdfLabel != "ep") {
                    continue;
                  }
                  std::string fullLabel = (pdfLabel != "") ? ("_" + pdfLabel) : "";
                  fillJetObservables(jetObservables.at("true" + fullLabel), std::get<1>(jetsTrue), jetR, kinematics);
                  fillJetObservables(jetObservables.at("true_charged" + fullLabel), std::get<1>(jetsTrueCharged), jetR, kinematics);
                  fillJetObservables(jetObservables.at("charged" + fullLabel), std::get<1>(jetsTrackRec), jetR, kinematics);
                  fillJetObservables(jetObservables.at("calo" + fullLabel), std::get<1>(jetsCaloRec), jetR, kinematics);
                  fillJetObservables(jetObservables.at("full" + fullLabel), std::get<1>(jetsFullRec), jetR, kinematics);
                }

                fillHadronObservables(jetObservables.at("true_ep"));

                if (jetR == jetRForJES) {
                  jetresolutionhistos(jetsTrackRec, jetsTrueCharged, 0, jetR);
                  jetresolutionhistos(jetsFullRec, jetsTrue, 1, jetR);
                  jetresolutionhistos(jetsHcalRec, jetsTrue, 2, jetR);
                  jetresolutionhistos(jetsCaloRec, jetsTrue, 3, jetR);
                  jetresolutionhistos(jetsAllRec, jetsTrue, 4, jetR);
                  jetresolutionhistos(jetsNoCluster,  jetsTrue,  5, jetR);
                  jetresolutionhistos(jetsEmcalRec,  jetsTrue,  6, jetR);
                }
              }
            }
            // TString jettype[njettypes] = {"track", "full","hcal","calo","all"};
        }
        if(tracksEnabled){
          if(verbosity>1) cout << "running trackingefficiency" << endl;
          trackingefficiency();
          if(verbosity>1) cout << "running trackingresolution" << endl;
          trackingresolution();
          if(verbosity>1) cout << "running trackingcomparison" << endl;
          trackingcomparison();
          if(verbosity>1) cout << "finished tracking studies" << endl;
        }
        if (runCaloRes){
          if(verbosity>1) cout << "running clusterstudies" << endl;
          clusterstudies();
          if(verbosity>1) std::cout << "running  caloresolutionhistos" << std::endl;
          caloresolutionhistos();
          if(verbosity>1) std::cout << "loop done ... next event" << std::endl;
        }
        if(tracksEnabled) trackmatchingstudies();

        clearClusterVectors();
        clearMCRecMatchVectors();

    } // event loop end
    if (_do_jetfinding) {
        std::cout << "saving event level observables\n";
        eventObservables.Write(outputDir.Data());
        std::cout << "saving jet observables\n";

        // The strategy here is that we want to overwrite previous outputs, but we otherwise want to group similar
        // results in the same file. To do this, we have each type with a separate output file write the first time
        // with "RECREATE", and then the rest write with "UPDATE"
        std::map<std::string, bool> firstWrite;
        for (const auto & pdfLabel : nPDFNames) {
          firstWrite.emplace(pdfLabel, true);
        }
        for (auto && [k, v] : jetObservables) {
          //std::cout << "Writing: " << k << "\n";
          std::string writeOption = "UPDATE";
          for (auto && [pdfLabel, first] : firstWrite) {
            if (k.find("_" + pdfLabel) != std::string::npos && first == true) {
              //std::cout << "First write for " << pdfLabel << "\n";
              std::cout << "Writing file for " << pdfLabel << "\n";
              writeOption = "RECREATE";
              first = false;
            }
          }
          v.Write(outputDir.Data(), writeOption);
        }
    }
    std::cout << "running jetresolutionhistosSave" << std::endl;
    jetresolutionhistosSave();
    if (runCaloRes){
      std::cout << "running caloresolutionhistosSave" << std::endl;
      caloresolutionhistossave();
      std::cout << "running clusterstudiesSave" << std::endl;
      clusterstudiesSave();
    }
    if(tracksEnabled){
      std::cout << "running trackingefficiencyhistosSave" << std::endl;
      trackingefficiencyhistosSave();
    }
    if(tracksEnabled) {
      std::cout << "running trackingresolutionhistosSave" << std::endl;
      trackingresolutionhistosSave();
    }
    if(tracksEnabled) {
      std::cout << "running trackingcomparisonhistosSave" << std::endl;
      trackingcomparisonhistosSave();
    }
    if(tracksEnabled){
      std::cout << "running hitstudiesSave" << std::endl;
      hitstudiesSave();
    }
    if(tracksEnabled) {
      std::cout << "running trackmatchingstudiesSave" << std::endl;
      trackmatchingstudiesSave();
    }
    std::cout << "running clusterizerSave" << std::endl;
    clusterizerSave();
    std::cout << "all done :)" << std::endl;
}
