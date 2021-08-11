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
        if(verbosity>1) std::cout << "event " << i << std::endl;

        // calculate useful quantities
        for(Int_t imc=0; imc<_nMCPart; imc++){
            TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
            _mcpart_Eta[imc]=truevec.Eta();
            int motherid = 0;
            for(Int_t imc2=0; imc2<_nMCPart; imc2++){
                if(_mcpart_ID[imc2]==_mcpart_ID_parent[imc]){
                    motherid=imc2;
                }
            }
            // if(_mcpart_PDG[imc]==111) cout << "pi0 found (" << _mcpart_ID[imc] << ") with mother: " << _mcpart_PDG[motherid] << " (" <<  _mcpart_ID_parent[imc] << ")" << endl;
        }
        float seed_E = 0.5;
        float aggregation_E = 0.1;
        // run clusterizers FHCAL
        if(do_reclus && _nTowers_FHCAL && caloEnabled[kFHCAL]){
            if(k3x3<_active_algo){
                if(verbosity>1) cout << "clusterizing 3x3 for FHCAL" << endl;
                runclusterizer(k3x3, kFHCAL,seed_E, aggregation_E, primaryTrackSource,
                    _nclusters_3x3_FHCAL,
                    _clusters_3x3_FHCAL_E,
                    _clusters_3x3_FHCAL_Eta,
                    _clusters_3x3_FHCAL_Phi,
                    _clusters_3x3_FHCAL_M02,
                    _clusters_3x3_FHCAL_M20,
                    _clusters_3x3_FHCAL_isMatched,
                    _clusters_3x3_FHCAL_NTower,
                    _clusters_3x3_FHCAL_trueID,
                    _clusters_3x3_FHCAL_NtrueID,
                    _clusters_3x3_FHCAL_X,
                    _clusters_3x3_FHCAL_Y,
                    _clusters_3x3_FHCAL_Z);
            }
            if(k5x5<_active_algo){
                if(verbosity>1) cout << "clusterizing 5x5 for FHCAL" << endl;
                runclusterizer(k5x5, kFHCAL,seed_E, aggregation_E, primaryTrackSource,
                    _nclusters_5x5_FHCAL,
                    _clusters_5x5_FHCAL_E,
                    _clusters_5x5_FHCAL_Eta,
                    _clusters_5x5_FHCAL_Phi,
                    _clusters_5x5_FHCAL_M02,
                    _clusters_5x5_FHCAL_M20,
                    _clusters_5x5_FHCAL_isMatched,
                    _clusters_5x5_FHCAL_NTower,
                    _clusters_5x5_FHCAL_trueID,
                    _clusters_5x5_FHCAL_NtrueID,
                    _clusters_5x5_FHCAL_X,
                    _clusters_5x5_FHCAL_Y,
                    _clusters_5x5_FHCAL_Z);
            }
            if(kV3<_active_algo){
                if(verbosity>1) cout << "clusterizing V3 for FHCAL" << endl;
                runclusterizer(kV3, kFHCAL,seed_E, aggregation_E, primaryTrackSource,
                    _nclusters_V3_FHCAL,
                    _clusters_V3_FHCAL_E,
                    _clusters_V3_FHCAL_Eta,
                    _clusters_V3_FHCAL_Phi,
                    _clusters_V3_FHCAL_M02,
                    _clusters_V3_FHCAL_M20,
                    _clusters_V3_FHCAL_isMatched,
                    _clusters_V3_FHCAL_NTower,
                    _clusters_V3_FHCAL_trueID,
                    _clusters_V3_FHCAL_NtrueID,
                    _clusters_V3_FHCAL_X,
                    _clusters_V3_FHCAL_Y,
                    _clusters_V3_FHCAL_Z);
            }
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for FHCAL" << endl;
                runclusterizer(kMA, kFHCAL,seed_E, aggregation_E, primaryTrackSource,
                    _nclusters_MA_FHCAL,
                    _clusters_MA_FHCAL_E,
                    _clusters_MA_FHCAL_Eta,
                    _clusters_MA_FHCAL_Phi,
                    _clusters_MA_FHCAL_M02,
                    _clusters_MA_FHCAL_M20,
                    _clusters_MA_FHCAL_isMatched,
                    _clusters_MA_FHCAL_NTower,
                    _clusters_MA_FHCAL_trueID,
                    _clusters_MA_FHCAL_NtrueID,
                    _clusters_MA_FHCAL_X,
                    _clusters_MA_FHCAL_Y,
                    _clusters_MA_FHCAL_Z);
            }
            if(kC3<_active_algo){
                if(verbosity>1) cout << "clusterizing C3 for FHCAL" << endl;
                runclusterizer(kC3, kFHCAL,seed_E, aggregation_E, primaryTrackSource,
                    _nclusters_C3_FHCAL,
                    _clusters_C3_FHCAL_E,
                    _clusters_C3_FHCAL_Eta,
                    _clusters_C3_FHCAL_Phi,
                    _clusters_C3_FHCAL_M02,
                    _clusters_C3_FHCAL_M20,
                    _clusters_C3_FHCAL_isMatched,
                    _clusters_C3_FHCAL_NTower,
                    _clusters_C3_FHCAL_trueID,
                    _clusters_C3_FHCAL_NtrueID,
                    _clusters_C3_FHCAL_X,
                    _clusters_C3_FHCAL_Y,
                    _clusters_C3_FHCAL_Z);
            }
            if(kC5<_active_algo){
                if(verbosity>1) cout << "clusterizing C5 for FHCAL" << endl;
                runclusterizer(kC5, kFHCAL,seed_E, aggregation_E, primaryTrackSource,
                    _nclusters_C5_FHCAL,
                    _clusters_C5_FHCAL_E,
                    _clusters_C5_FHCAL_Eta,
                    _clusters_C5_FHCAL_Phi,
                    _clusters_C5_FHCAL_M02,
                    _clusters_C5_FHCAL_M20,
                    _clusters_C5_FHCAL_isMatched,
                    _clusters_C5_FHCAL_NTower,
                    _clusters_C5_FHCAL_trueID,
                    _clusters_C5_FHCAL_NtrueID,
                    _clusters_C5_FHCAL_X,
                    _clusters_C5_FHCAL_Y,
                    _clusters_C5_FHCAL_Z);
            }
        } else {
            _nclusters_3x3_FHCAL = 0;
            _nclusters_5x5_FHCAL = 0;
            _nclusters_V3_FHCAL = 0;
            _nclusters_MA_FHCAL = 0;
            _nclusters_C3_FHCAL = 0;
            _nclusters_C5_FHCAL = 0;
        }
        // run clusterizers FEMC
        if(do_reclus && _nTowers_FEMC && caloEnabled[kFEMC]){
          Float_t seed_E_FEMC         = 0.1;
          Float_t aggregation_E_FEMC  = 0.005;
            if(k3x3<_active_algo){
                if(verbosity>1) cout << "clusterizing 3x3 for FEMC" << endl;
                runclusterizer(k3x3, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource,
                    _nclusters_3x3_FEMC,
                    _clusters_3x3_FEMC_E,
                    _clusters_3x3_FEMC_Eta,
                    _clusters_3x3_FEMC_Phi,
                    _clusters_3x3_FEMC_M02,
                    _clusters_3x3_FEMC_M20,
                    _clusters_3x3_FEMC_isMatched,
                    _clusters_3x3_FEMC_NTower,
                    _clusters_3x3_FEMC_trueID,
                    _clusters_3x3_FEMC_NtrueID,
                    _clusters_3x3_FEMC_X,
                    _clusters_3x3_FEMC_Y,
                    _clusters_3x3_FEMC_Z);
            }
            if(k5x5<_active_algo){
                if(verbosity>1) cout << "clusterizing 5x5 for FEMC" << endl;
                runclusterizer(k5x5, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource,
                    _nclusters_5x5_FEMC,
                    _clusters_5x5_FEMC_E,
                    _clusters_5x5_FEMC_Eta,
                    _clusters_5x5_FEMC_Phi,
                    _clusters_5x5_FEMC_M02,
                    _clusters_5x5_FEMC_M20,
                    _clusters_5x5_FEMC_isMatched,
                    _clusters_5x5_FEMC_NTower,
                    _clusters_5x5_FEMC_trueID,
                    _clusters_5x5_FEMC_NtrueID,
                    _clusters_5x5_FEMC_X,
                    _clusters_5x5_FEMC_Y,
                    _clusters_5x5_FEMC_Z);
            }
            if(kV3<_active_algo){
                if(verbosity>1) cout << "clusterizing V3 for FEMC" << endl;
                runclusterizer(kV3, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource,
                    _nclusters_V3_FEMC,
                    _clusters_V3_FEMC_E,
                    _clusters_V3_FEMC_Eta,
                    _clusters_V3_FEMC_Phi,
                    _clusters_V3_FEMC_M02,
                    _clusters_V3_FEMC_M20,
                    _clusters_V3_FEMC_isMatched,
                    _clusters_V3_FEMC_NTower,
                    _clusters_V3_FEMC_trueID,
                    _clusters_V3_FEMC_NtrueID,
                    _clusters_V3_FEMC_X,
                    _clusters_V3_FEMC_Y,
                    _clusters_V3_FEMC_Z);
            }
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for FEMC" << endl;
                runclusterizer(kMA, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource,
                    _nclusters_MA_FEMC,
                    _clusters_MA_FEMC_E,
                    _clusters_MA_FEMC_Eta,
                    _clusters_MA_FEMC_Phi,
                    _clusters_MA_FEMC_M02,
                    _clusters_MA_FEMC_M20,
                    _clusters_MA_FEMC_isMatched,
                    _clusters_MA_FEMC_NTower,
                    _clusters_MA_FEMC_trueID,
                    _clusters_MA_FEMC_NtrueID,
                    _clusters_MA_FEMC_X,
                    _clusters_MA_FEMC_Y,
                    _clusters_MA_FEMC_Z);
            }
            if(kC3<_active_algo){
                if(verbosity>1) cout << "clusterizing C3 for FEMC" << endl;
                runclusterizer(kC3, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource,
                    _nclusters_C3_FEMC,
                    _clusters_C3_FEMC_E,
                    _clusters_C3_FEMC_Eta,
                    _clusters_C3_FEMC_Phi,
                    _clusters_C3_FEMC_M02,
                    _clusters_C3_FEMC_M20,
                    _clusters_C3_FEMC_isMatched,
                    _clusters_C3_FEMC_NTower,
                    _clusters_C3_FEMC_trueID,
                    _clusters_C3_FEMC_NtrueID,
                    _clusters_C3_FEMC_X,
                    _clusters_C3_FEMC_Y,
                    _clusters_C3_FEMC_Z);
            }
            if(kC5<_active_algo){
                if(verbosity>1) cout << "clusterizing C5 for FEMC" << endl;
                runclusterizer(kC5, kFEMC,seed_E_FEMC, aggregation_E_FEMC, primaryTrackSource,
                    _nclusters_C5_FEMC,
                    _clusters_C5_FEMC_E,
                    _clusters_C5_FEMC_Eta,
                    _clusters_C5_FEMC_Phi,
                    _clusters_C5_FEMC_M02,
                    _clusters_C5_FEMC_M20,
                    _clusters_C5_FEMC_isMatched,
                    _clusters_C5_FEMC_NTower,
                    _clusters_C5_FEMC_trueID,
                    _clusters_C5_FEMC_NtrueID,
                    _clusters_C5_FEMC_X,
                    _clusters_C5_FEMC_Y,
                    _clusters_C5_FEMC_Z);
            }
        } else {
            _nclusters_3x3_FEMC = 0;
            _nclusters_5x5_FEMC = 0;
            _nclusters_V3_FEMC = 0;
            _nclusters_MA_FEMC = 0;
            _nclusters_C3_FEMC = 0;
            _nclusters_C5_FEMC = 0;
        }

        if(do_reclus && _nTowers_CEMC && caloEnabled[kCEMC]){
            float seed_E_CEMC = 0.5;
            float aggregation_E_CEMC = 0.1;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for CEMC" << endl;
                runclusterizer(kMA, kCEMC,seed_E_CEMC, aggregation_E_CEMC, primaryTrackSource,
                    _nclusters_MA_CEMC,
                    _clusters_MA_CEMC_E,
                    _clusters_MA_CEMC_Eta,
                    _clusters_MA_CEMC_Phi,
                    _clusters_MA_CEMC_M02,
                    _clusters_MA_CEMC_M20,
                    _clusters_MA_CEMC_isMatched,
                    _clusters_MA_CEMC_NTower,
                    _clusters_MA_CEMC_trueID,
                    _clusters_MA_CEMC_NtrueID,
                    _clusters_MA_CEMC_X,
                    _clusters_MA_CEMC_Y,
                    _clusters_MA_CEMC_Z);
            }
        } else _nclusters_MA_CEMC = 0;

        if(do_reclus && _nTowers_HCALIN && caloEnabled[kHCALIN]){
            float seed_E_HCALIN = 0.2;
            float aggregation_E_HCALIN = 0.05;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for HCALIN" << endl;
                runclusterizer(kMA, kHCALIN,seed_E_HCALIN, aggregation_E_HCALIN, primaryTrackSource,
                    _nclusters_MA_HCALIN,
                    _clusters_MA_HCALIN_E,
                    _clusters_MA_HCALIN_Eta,
                    _clusters_MA_HCALIN_Phi,
                    _clusters_MA_HCALIN_M02,
                    _clusters_MA_HCALIN_M20,
                    _clusters_MA_HCALIN_isMatched,
                    _clusters_MA_HCALIN_NTower,
                    _clusters_MA_HCALIN_trueID,
                    _clusters_MA_HCALIN_NtrueID,
                    _clusters_MA_HCALIN_X,
                    _clusters_MA_HCALIN_Y,
                    _clusters_MA_HCALIN_Z);
            }
        } else _nclusters_MA_HCALIN = 0;

        if(do_reclus && _nTowers_HCALOUT && caloEnabled[kHCALOUT]){
            float seed_E_HCALOUT = 0.5;
            float aggregation_E_HCALOUT = 0.1;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for HCALOUT" << endl;
                runclusterizer(kMA, kHCALOUT,seed_E_HCALOUT, aggregation_E_HCALOUT, primaryTrackSource,
                    _nclusters_MA_HCALOUT,
                    _clusters_MA_HCALOUT_E,
                    _clusters_MA_HCALOUT_Eta,
                    _clusters_MA_HCALOUT_Phi,
                    _clusters_MA_HCALOUT_M02,
                    _clusters_MA_HCALOUT_M20,
                    _clusters_MA_HCALOUT_isMatched,
                    _clusters_MA_HCALOUT_NTower,
                    _clusters_MA_HCALOUT_trueID,
                    _clusters_MA_HCALOUT_NtrueID,
                    _clusters_MA_HCALOUT_X,
                    _clusters_MA_HCALOUT_Y,
                    _clusters_MA_HCALOUT_Z);
            }
        } else _nclusters_MA_HCALOUT = 0;

        if(do_reclus && _nTowers_EEMCG && caloEnabled[kEEMCG]){
            float seed_E_EEMCG = 0.1;
            float aggregation_E_EEMCG = 0.01;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for EEMCG" << endl;
                runclusterizer(kMA, kEEMCG,seed_E_EEMCG, aggregation_E_EEMCG, primaryTrackSource,
                    _nclusters_MA_EEMCG,
                    _clusters_MA_EEMCG_E,
                    _clusters_MA_EEMCG_Eta,
                    _clusters_MA_EEMCG_Phi,
                    _clusters_MA_EEMCG_M02,
                    _clusters_MA_EEMCG_M20,
                    _clusters_MA_EEMCG_isMatched,
                    _clusters_MA_EEMCG_NTower,
                    _clusters_MA_EEMCG_trueID,
                    _clusters_MA_EEMCG_NtrueID,
                    _clusters_MA_EEMCG_X,
                    _clusters_MA_EEMCG_Y,
                    _clusters_MA_EEMCG_Z);
            }
        } else _nclusters_MA_EEMCG = 0;

        if(do_reclus && _nTowers_BECAL  && caloEnabled[kBECAL]){
            float seed_E_BECAL = 0.1;
            float aggregation_E_BECAL = 0.01;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for BECAL" << endl;
                runclusterizer(kMA, kBECAL,seed_E_BECAL, aggregation_E_BECAL, primaryTrackSource,
                    _nclusters_MA_BECAL,
                    _clusters_MA_BECAL_E,
                    _clusters_MA_BECAL_Eta,
                    _clusters_MA_BECAL_Phi,
                    _clusters_MA_BECAL_M02,
                    _clusters_MA_BECAL_M20,
                    _clusters_MA_BECAL_isMatched,
                    _clusters_MA_BECAL_NTower,
                    _clusters_MA_BECAL_trueID,
                    _clusters_MA_BECAL_NtrueID,
                    _clusters_MA_BECAL_X,
                    _clusters_MA_BECAL_Y,
                    _clusters_MA_BECAL_Z);
            }
        } else _nclusters_MA_BECAL = 0;

        if(do_reclus && _nTowers_EEMC  && caloEnabled[kEEMC]){
            float seed_E_EEMC = 0.1;
            float aggregation_E_EEMC = 0.05;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for EEMC" << endl;
                runclusterizer(kMA, kEEMC,seed_E_EEMC, aggregation_E_EEMC, primaryTrackSource,
                    _nclusters_MA_EEMC,
                    _clusters_MA_EEMC_E,
                    _clusters_MA_EEMC_Eta,
                    _clusters_MA_EEMC_Phi,
                    _clusters_MA_EEMC_M02,
                    _clusters_MA_EEMC_M20,
                    _clusters_MA_EEMC_isMatched,
                    _clusters_MA_EEMC_NTower,
                    _clusters_MA_EEMC_trueID,
                    _clusters_MA_EEMC_NtrueID,
                    _clusters_MA_EEMC_X,
                    _clusters_MA_EEMC_Y,
                    _clusters_MA_EEMC_Z);
            }
        } else _nclusters_MA_EEMC = 0;

        if(do_reclus && _nTowers_LFHCAL && caloEnabled[kLFHCAL]){
            float seed_E_LFHCAL = 0.1;
            float aggregation_E_LFHCAL = 0.001;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for LFHCAL" << endl;
                runclusterizer(kMA, kLFHCAL,seed_E_LFHCAL, aggregation_E_LFHCAL, primaryTrackSource,
                    _nclusters_MA_LFHCAL,
                    _clusters_MA_LFHCAL_E,
                    _clusters_MA_LFHCAL_Eta,
                    _clusters_MA_LFHCAL_Phi,
                    _clusters_MA_LFHCAL_M02,
                    _clusters_MA_LFHCAL_M20,
                    _clusters_MA_LFHCAL_isMatched,
                    _clusters_MA_LFHCAL_NTower,
                    _clusters_MA_LFHCAL_trueID,
                    _clusters_MA_LFHCAL_NtrueID,
                    _clusters_MA_LFHCAL_X,
                    _clusters_MA_LFHCAL_Y,
                    _clusters_MA_LFHCAL_Z);
            }
        } else _nclusters_MA_LFHCAL = 0;

        if(do_reclus && _nTowers_EHCAL && caloEnabled[kEHCAL]){
            float seed_E_EHCAL = 0.01;
            float aggregation_E_EHCAL = 0.005;
            if(kMA<_active_algo){
                if(verbosity>1) cout << "clusterizing MA for EHCAL" << endl;
                runclusterizer(kMA, kEHCAL,seed_E_EHCAL, aggregation_E_EHCAL, primaryTrackSource,
                    _nclusters_MA_EHCAL,
                    _clusters_MA_EHCAL_E,
                    _clusters_MA_EHCAL_Eta,
                    _clusters_MA_EHCAL_Phi,
                    _clusters_MA_EHCAL_M02,
                    _clusters_MA_EHCAL_M20,
                    _clusters_MA_EHCAL_isMatched,
                    _clusters_MA_EHCAL_NTower,
                    _clusters_MA_EHCAL_trueID,
                    _clusters_MA_EHCAL_NtrueID,
                    _clusters_MA_EHCAL_X,
                    _clusters_MA_EHCAL_Y,
                    _clusters_MA_EHCAL_Z);
            }
        } else _nclusters_MA_EHCAL = 0;

        if(do_reclus && _nTowers_DRCALO && caloEnabled[kDRCALO]){ //do_V1clusterizerDRCALO
            float seed_E_DRCALO = 0.3;
            float aggregation_E_DRCALO = 0.05;
            if(verbosity>1) cout << "clusterizing V1 for DRCALO" << endl;
            runclusterizer(kV1, kDRCALO,seed_E_DRCALO, aggregation_E_DRCALO, primaryTrackSource,
                _nclusters_V1_DRCALO,
                _clusters_V1_DRCALO_E,
                _clusters_V1_DRCALO_Eta,
                _clusters_V1_DRCALO_Phi,
                _clusters_V1_DRCALO_M02,
                _clusters_V1_DRCALO_M20,
                _clusters_V1_DRCALO_isMatched,
                _clusters_V1_DRCALO_NTower,
                _clusters_V1_DRCALO_trueID,
                _clusters_V1_DRCALO_NtrueID,
                _clusters_V1_DRCALO_X,
                _clusters_V1_DRCALO_Y,
                _clusters_V1_DRCALO_Z);
        } else _nclusters_V1_DRCALO = 0;
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

        // ANCHOR Track projections loop variables:
        // float* _track_ProjTrackID[iproj]
        // int* _track_ProjLayer[iproj]
        // float* _track_Proj_x[iproj]
        // float* _track_Proj_y[iproj]
        // float* _track_Proj_z[iproj]
        // float* _track_Proj_t[iproj]
        // float* _track_Proj_true_x[iproj]
        // float* _track_Proj_true_y[iproj]
        // float* _track_Proj_true_z[iproj]
        // float* _track_Proj_true_t[iproj]
        // for(Int_t iproj=0; iproj<_nProjections; iproj++){
        //     if(verbosity>1) std::cout << "\tProjection: proj " << iproj << "\tin layer " << _track_ProjLayer[iproj] << "\twith X = " << _track_Proj_x[iproj] << " cm" << std::endl;
        // }

        // ANCHOR FHCAL cluster loop variables:
        // float* _clusters_FHCAL_E[iclus];
        // float* _clusters_FHCAL_Eta[iclus];
        // float* _clusters_FHCAL_Phi[iclus];
        // int* _clusters_FHCAL_NTower[iclus];
        // float* _clusters_FHCAL_trueID[iclus];
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
        // for(Int_t iclus=0; iclus<_nclusters_FHCAL; iclus++){
        //     if(verbosity>1) std::cout << "\tcls " << iclus << "\tE " << _clusters_FHCAL_E[iclus] << "\tEta " << _clusters_FHCAL_Eta[iclus] << "\tPhi " << _clusters_FHCAL_Phi[iclus] << "\tntowers: " << _clusters_FHCAL_NTower[iclus] << "\ttrueID: " << _clusters_FHCAL_trueID[iclus] << std::endl;
        // }
        // apply calibration if desired
        if(kV1<_active_algo){
            // TODO: Should kV1 in the calibration actually be _active_algo?
            for(Int_t iclus=0; iclus<_nclusters_FHCAL; iclus++){
                _clusters_FHCAL_trueID[iclus] = GetCorrectMCArrayEntry(_clusters_FHCAL_trueID[iclus]);
                if(_doClusterECalibration){
                    _clusters_FHCAL_E[iclus]/=getCalibrationValue(_clusters_FHCAL_E[iclus], kFHCAL, kV1);
                }
            }
            for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
                _clusters_FEMC_trueID[iclus] = GetCorrectMCArrayEntry(_clusters_FEMC_trueID[iclus]);
                if(_doClusterECalibration){
                    _clusters_FEMC_E[iclus]/=getCalibrationValue(_clusters_FEMC_E[iclus], kFEMC, kV1);
                }
            }
        }
        if(do_reclus && kMA<_active_algo && _do_jetfinding){
            fillHCalClustersIntoJetFindingInputs(_nclusters_MA_FHCAL, _clusters_MA_FHCAL_E, _clusters_MA_FHCAL_Eta, _clusters_MA_FHCAL_Phi, _clusters_MA_FHCAL_isMatched,
                                                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                                                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                                                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz);
        }

        // ANCHOR FEMC cluster loop variables:
        // float* _clusters_FEMC_E[iclus];
        // float* _clusters_FEMC_Eta[iclus];
        // float* _clusters_FEMC_Phi[iclus];
        // int* _clusters_FEMC_NTower[iclus];
        // float* _clusters_FEMC_trueID[iclus];
        // for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
        //     if(verbosity>1) std::cout << "\tFEMC:  cluster " << iclus << "\twith E = " << _clusters_FEMC_E[iclus] << " GeV" << std::endl;
        // }
        if(do_reclus && kMA<_active_algo && _do_jetfinding){
            fillECalClustersIntoJetFindingInputs(_nclusters_MA_FEMC, _clusters_MA_FEMC_E, _clusters_MA_FEMC_Eta, _clusters_MA_FEMC_Phi, _clusters_MA_FEMC_isMatched,
                                                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                                                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                                                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                                                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz);
        }

        // Add rest of calorimeter clusters to jet finder inputs
        if(do_reclus && kMA<_active_algo && _do_jetfinding) {
            // EEMC
            fillECalClustersIntoJetFindingInputs(
                _nclusters_MA_EEMCG, _clusters_MA_EEMCG_E, _clusters_MA_EEMCG_Eta, _clusters_MA_EEMCG_Phi, _clusters_MA_EEMCG_isMatched,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // EEMCG
            fillECalClustersIntoJetFindingInputs(
                _nclusters_MA_EEMCG, _clusters_MA_EEMCG_E, _clusters_MA_EEMCG_Eta, _clusters_MA_EEMCG_Phi, _clusters_MA_EEMCG_isMatched,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // EHCAL
            fillHCalClustersIntoJetFindingInputs(
                _nclusters_MA_EHCAL, _clusters_MA_EHCAL_E, _clusters_MA_EHCAL_Eta, _clusters_MA_EHCAL_Phi, _clusters_MA_EHCAL_isMatched,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
            // CEMC
            fillECalClustersIntoJetFindingInputs(
                _nclusters_MA_CEMC, _clusters_MA_CEMC_E, _clusters_MA_CEMC_Eta, _clusters_MA_CEMC_Phi, _clusters_MA_CEMC_isMatched,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // BECAL
            fillECalClustersIntoJetFindingInputs(
                _nclusters_MA_BECAL, _clusters_MA_BECAL_E, _clusters_MA_BECAL_Eta, _clusters_MA_BECAL_Phi, _clusters_MA_BECAL_isMatched,
                jetf_emcal_E, jetf_emcal_px, jetf_emcal_py, jetf_emcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz,
                jetf_full_E, jetf_full_px, jetf_full_py, jetf_full_pz
            );
            // HCALIN
            fillHCalClustersIntoJetFindingInputs(
                _nclusters_MA_HCALIN, _clusters_MA_HCALIN_E, _clusters_MA_HCALIN_Eta, _clusters_MA_HCALIN_Phi, _clusters_MA_HCALIN_isMatched,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
            // HCALOUT
            fillHCalClustersIntoJetFindingInputs(
                _nclusters_MA_HCALOUT, _clusters_MA_HCALOUT_E, _clusters_MA_HCALOUT_Eta, _clusters_MA_HCALOUT_Phi, _clusters_MA_HCALOUT_isMatched,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
            // LFHCAL
            fillHCalClustersIntoJetFindingInputs(
                _nclusters_MA_LFHCAL, _clusters_MA_LFHCAL_E, _clusters_MA_LFHCAL_Eta, _clusters_MA_LFHCAL_Phi, _clusters_MA_LFHCAL_isMatched,
                jetf_hcal_E, jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz,
                jetf_calo_E, jetf_calo_px, jetf_calo_py, jetf_calo_pz,
                jetf_all_E, jetf_all_px, jetf_all_py, jetf_all_pz
            );
        }

        // ANCHOR FHCAL tower loop variables:
        // float* _tower_FHCAL_E[itwrH]
        // int* _tower_FHCAL_iEta[itwrH]
        // int* _tower_FHCAL_iPhi[itwrH]
        // float* _tower_FHCAL_trueID[itwrH]
        // for(Int_t itwrH=0; itwrH<_nTowers_FHCAL; itwrH++){
        //     if(verbosity>1) std::cout << "\tFHCAL: tower " << itwrH << "\twith E = " << _tower_FHCAL_E[itwrH] << " GeV" << std::endl;
        // }

        // ANCHOR FEMC cluster loop variables:
        // float* _tower_FEMC_E[itwrH]
        // int* _tower_FEMC_iEta[itwrH]
        // int* _tower_FEMC_iPhi[itwrH]
        // float* _tower_FEMC_trueID[itwrH]
        // for(Int_t itwrE=0; itwrE<_nTowers_FEMC; itwrE++){
        //     if(verbosity>1) std::cout << "\tFEMC:  tower " << itwrE << "\twith E = " << _tower_FEMC_E[itwrE] << " GeV" << std::endl;
        // }

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
                for (uint32_t i = 0; i < _nTowers_FEMC; i++) {
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
