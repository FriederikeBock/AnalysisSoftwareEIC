#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include "treeProcessing.h"
#include "event_utils.cxx"
#include "jet_finder.cxx"
#include "jet_observables.cxx"
#include "caloheader.h"
#include "clusterizer.cxx"

#include "jetresolutionhistos.cxx"
#include "resolutionhistos.cxx"
#include "clusterstudies.cxx"
#include "trackingefficiency.cxx"
#include "hitstudies.cxx"
#include "trackmatchingstudies.cxx"

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TVector3.h>


#include <iostream>
#include <fstream>

void treeProcessing(
    TString inFile              = "",
    TString inFileGeometry      = "geometry.root",
    TString addOutputName       = "",
    bool do_reclus              = true,
    bool do_jetfinding          = false,
    // Double_t maxNEvent = 1e5,
    bool hasTiming              = true,
    bool isALLSILICON           = true,
    Double_t maxNEvent          = -1,
    Int_t verbosity             = 0,
    bool doCalibration          = false,
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0,
    std::string jetAlgorithm    = "anti-kt",
    double jetR                 = 0.5,
    double tracked_jet_max_pT   = 30
){
    // make output directory
    TString dateForOutput = ReturnDateStr();
    outputDir = Form("treeProcessing/%s",addOutputName.Data());
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
        std::cout << "loading a list of files" << std::endl;
        std::ifstream files(inFile);
        std::string filePath;

        while (std::getline(files, filePath)) {
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

    for (Int_t c = 0; c < 11; c++){
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
    auto eventObservables = EventObservables();
    auto jetObservablesTrue = JetObservables{JetType_t::full, "true"};
    auto jetObservablesTrueCharged = JetObservables{JetType_t::charged, "true"};
    auto jetObservablesCharged = JetObservables{JetType_t::charged};
    auto jetObservablesCalo = JetObservables{JetType_t::calo};
    auto jetObservablesFull = JetObservables{JetType_t::full};

    if (_do_jetfinding) {
        // Event level
        eventObservables.Init();

        // TODO: Add fully to arguments
        std::vector<double> jetRParameters = {0.5};
        jetObservablesTrue.Init(jetRParameters);
        jetObservablesTrueCharged.Init(jetRParameters);
        jetObservablesCharged.Init(jetRParameters);
        jetObservablesCalo.Init(jetRParameters);
        jetObservablesFull.Init(jetRParameters);
    }

    _nEventsTree=0;
    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {
        // load current event
        tt_event->GetEntry(i);
        _nEventsTree++;

        // processing progress info
        if(i>0 && nEntriesTree>100 && i%(nEntriesTree/(20))==0) std::cout << "//processed " << 100*(i)/nEntriesTree << "%"  << std::endl;
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

        
        float seed_E = 0.5;
        float aggregation_E = 0.1;
        // run clusterizers FHCAL
        if(do_reclus && _nTowers_FHCAL && caloEnabled[kFHCAL]){
            if(k3x3<_active_algo){
                if(verbosity>1) cout << "clusterizing 3x3 for FHCAL" << endl;
                runclusterizer(k3x3, kFHCAL,seed_E, aggregation_E, primaryTrackSource);
            }
            if(k5x5<_active_algo){
                if(verbosity>1) cout << "clusterizing 5x5 for FHCAL" << endl;
                runclusterizer(k5x5, kFHCAL,seed_E, aggregation_E, primaryTrackSource);
            }
            if(kV3<_active_algo){
                if(verbosity>1) cout << "clusterizing V3 for FHCAL" << endl;
                runclusterizer(kV3, kFHCAL,seed_E, aggregation_E, primaryTrackSource);
            }
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for FHCAL" << endl;
                runclusterizer(kMA, kFHCAL,seed_E, aggregation_E, primaryTrackSource);
            }
            if(kC3<_active_algo){
                if(verbosity>1) cout << "clusterizing C3 for FHCAL" << endl;
                runclusterizer(kC3, kFHCAL,seed_E, aggregation_E, primaryTrackSource);
            }
            if(kC5<_active_algo){
                if(verbosity>1) cout << "clusterizing C5 for FHCAL" << endl;
                runclusterizer(kC5, kFHCAL,seed_E, aggregation_E, primaryTrackSource);
            }
        } 
        
        // run clusterizers FEMC
        if(do_reclus && _nTowers_FEMC && caloEnabled[kFEMC]){
          Float_t seed_E_FEMC         = 0.1;
          Float_t aggregation_E_FEMC  = 0.005;
            if(k3x3<_active_algo){
                if(verbosity>1) cout << "clusterizing 3x3 for FEMC" << endl;
                runclusterizer(k3x3, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource);
            }
            if(k5x5<_active_algo){
                if(verbosity>1) cout << "clusterizing 5x5 for FEMC" << endl;
                runclusterizer(k5x5, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource);
            }
            if(kV3<_active_algo){
                if(verbosity>1) cout << "clusterizing V3 for FEMC" << endl;
                runclusterizer(kV3, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource);
            }
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for FEMC" << endl;
                runclusterizer(kMA, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource);
            }
            if(kC3<_active_algo){
                if(verbosity>1) cout << "clusterizing C3 for FEMC" << endl;
                runclusterizer(kC3, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource);
            }
            if(kC5<_active_algo){
                if(verbosity>1) cout << "clusterizing C5 for FEMC" << endl;
                runclusterizer(kC5, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_CEMC && caloEnabled[kCEMC]){
            float seed_E_CEMC = 0.5;
            float aggregation_E_CEMC = 0.1;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for CEMC" << endl;
                runclusterizer(kMA, kCEMC,seed_E_CEMC, aggregation_E_CEMC, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_HCALIN && caloEnabled[kHCALIN]){
            float seed_E_HCALIN = 0.2;
            float aggregation_E_HCALIN = 0.05;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for HCALIN" << endl;
                runclusterizer(kMA, kHCALIN,seed_E_HCALIN, aggregation_E_HCALIN, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_HCALOUT && caloEnabled[kHCALOUT]){
            float seed_E_HCALOUT = 0.5;
            float aggregation_E_HCALOUT = 0.1;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for HCALOUT" << endl;
                runclusterizer(kMA, kHCALOUT,seed_E_HCALOUT, aggregation_E_HCALOUT, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_EEMCG && caloEnabled[kEEMCG]){
            float seed_E_EEMCG = 0.1;
            float aggregation_E_EEMCG = 0.01;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for EEMCG" << endl;
                runclusterizer(kMA, kEEMCG,seed_E_EEMCG, aggregation_E_EEMCG, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_BECAL  && caloEnabled[kBECAL]){
            float seed_E_BECAL = 0.1;
            float aggregation_E_BECAL = 0.01;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for BECAL" << endl;
                runclusterizer(kMA, kBECAL,seed_E_BECAL, aggregation_E_BECAL, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_EEMC  && caloEnabled[kEEMC]){
            float seed_E_EEMC = 0.1;
            float aggregation_E_EEMC = 0.05;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for EEMC" << endl;
                runclusterizer(kMA, kEEMC,seed_E_EEMC, aggregation_E_EEMC, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_LFHCAL && caloEnabled[kLFHCAL]){
            float seed_E_LFHCAL = 0.1;
            float aggregation_E_LFHCAL = 0.001;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for LFHCAL" << endl;
                runclusterizer(kMA, kLFHCAL,seed_E_LFHCAL, aggregation_E_LFHCAL, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_EHCAL && caloEnabled[kEHCAL]){
            float seed_E_EHCAL = 0.01;
            float aggregation_E_EHCAL = 0.005;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for EHCAL" << endl;
                runclusterizer(kMA, kEHCAL,seed_E_EHCAL, aggregation_E_EHCAL, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_DRCALO && caloEnabled[kDRCALO]){ //do_V1clusterizerDRCALO
            float seed_E_DRCALO = 0.3;
            float aggregation_E_DRCALO = 0.05;
            if(verbosity>1) cout << "clusterizing V1 for DRCALO" << endl;
            runclusterizer(kV1, kDRCALO,seed_E_DRCALO, aggregation_E_DRCALO, primaryTrackSource);
        } 
        if((do_reclus) && verbosity>1) cout << "done with clusterization!" << endl;

        // ANCHOR Hits loop variables:
        // float* _hits_x[ihit]
        // float* _hits_y[ihit]
        // float* _hits_z[ihit]
        // float* _hits_t[ihit]
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
            fillHCalClustersIntoJetFindingInputs( fwdCaloID, kMA,
                                                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                                                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                                                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz);
        }

        // ANCHOR FEMC cluster loop variables:
        // for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
        //     if(verbosity>1) std::cout << "\tFEMC:  cluster " << iclus << "\twith E = " << (_clusters_calo[kV1][kFEMC].at(iclus)).cluster_E << " GeV" << std::endl;
        // }
        if(do_reclus && kMA<_active_algo && _do_jetfinding){
            fillECalClustersIntoJetFindingInputs( kFEMC, kMA,
                                                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                                                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                                                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                                                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz);
        }

        // Add rest of calorimeter clusters to jet finder inputs
        if(do_reclus && kMA<_active_algo && _do_jetfinding) {
            // EEMC
            fillECalClustersIntoJetFindingInputs(
                kEEMC, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // EEMCG
            fillECalClustersIntoJetFindingInputs(
                kEEMCG, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // EHCAL
            fillHCalClustersIntoJetFindingInputs(
                kEHCAL, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
            // CEMC
            fillECalClustersIntoJetFindingInputs(
                kCEMC, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // BECAL
            fillECalClustersIntoJetFindingInputs(
                kBECAL, kMA,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // HCALIN
            fillHCalClustersIntoJetFindingInputs(
                kHCALIN, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
            // HCALOUT
            fillHCalClustersIntoJetFindingInputs(
                kHCALOUT, kMA,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
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

            std::vector<double> nocluster_px;
            std::vector<double> nocluster_py;
            std::vector<double> nocluster_pz;
            std::vector<double> nocluster_E;
            if (true) {
                for (uint32_t i = 0; i < (uint32_t)_nTowers_FEMC; i++) {
                    if (_tower_FEMC_E[i] < seed_E) {
                        continue;
                    }
                    TVector3 twrPos = TowerPositionVectorFromIndicesGeometry(_tower_FEMC_iEta[i], _tower_FEMC_iPhi[i], kFEMC);
                    double pt = _tower_FEMC_E[i] / cosh(twrPos.Eta());
                    nocluster_px.push_back(pt * cos(twrPos.Phi()));
                    nocluster_py.push_back(pt * sin(twrPos.Phi()));
                    nocluster_pz.push_back(pt * sinh(twrPos.Eta()));
                    nocluster_E.push_back(_tower_FEMC_E[i]);
                }
            }
            // ANCHOR JET FINDING
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
            fillEventObservables(eventObservables, primaryTrackSource);
            fillJetObservables(jetObservablesTrue, std::get<1>(jetsTrue), jetR);
            fillJetObservables(jetObservablesTrueCharged, std::get<1>(jetsTrueCharged), jetR);
            fillJetObservables(jetObservablesCharged, std::get<1>(jetsTrackRec), jetR);
            fillJetObservables(jetObservablesCalo, std::get<1>(jetsCaloRec), jetR);
            fillJetObservables(jetObservablesFull, std::get<1>(jetsFullRec), jetR);

            fillHadronObservables(jetObservablesTrue);

            jetresolutionhistos(jetsTrackRec, jetsTrueCharged, 0, jetR);
            jetresolutionhistos(jetsFullRec, jetsTrue, 1, jetR);
            jetresolutionhistos(jetsHcalRec, jetsTrue, 2, jetR);
            jetresolutionhistos(jetsCaloRec, jetsTrue, 3, jetR);
            jetresolutionhistos(jetsAllRec, jetsTrue, 4, jetR);
            jetresolutionhistos(jetsNoCluster,  jetsTrue,  5, jetR);
            jetresolutionhistos(jetsEmcalRec,  jetsTrue,  6, jetR);
            // TString jettype[njettypes] = {"track", "full","hcal","calo","all"};
        }
        if(tracksEnabled){        
          if(verbosity>1) cout << "running trackingefficiency" << endl;
          trackingefficiency();
          if(verbosity>1) cout << "running trackingresolution" << endl;
          trackingresolution();
          if(verbosity>1) cout << "running trackingcomparison" << endl;
          // trackingcomparison();
          if(verbosity>1) cout << "finished tracking studies" << endl;
        }
        if(verbosity>1) cout << "running clusterstudies" << endl;
        clusterstudies();
        if(verbosity>1) std::cout << "running resolutionhistos" << std::endl;
        resolutionhistos();
        if(verbosity>1) std::cout << "loop done ... next event" << std::endl;
        if(tracksEnabled) trackmatchingstudies();
        
        clearClusterVectors();

    } // event loop end
    if (_do_jetfinding) {
        std::cout << "saving event level observables\n";
        eventObservables.Write(outputDir.Data());
        std::cout << "saving jet observables\n";
        jetObservablesTrue.Write(outputDir.Data(), "RECREATE");
        jetObservablesTrueCharged.Write(outputDir.Data());
        jetObservablesCharged.Write(outputDir.Data());
        jetObservablesCalo.Write(outputDir.Data());
        jetObservablesFull.Write(outputDir.Data());
    }
    std::cout << "running jetresolutionhistosSave" << std::endl;
    jetresolutionhistosSave();
    std::cout << "running resolutionhistosSave" << std::endl;
    resolutionhistosSave();
    std::cout << "running clusterstudiesSave" << std::endl;
    clusterstudiesSave();
    
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
