#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include "treeProcessing.h"
#include "event_utils.cxx"
#include "caloheader.h"
#include "clusterizer.cxx"

//#include "disreconstruction.cxx"
#include "trackmatchingstudies.cxx"
#include "electronpid.cxx"

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TVector3.h>


#include <iostream>
#include <fstream>

void epidTreeProcessing(
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
    bool doCalibration          = true,
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

    make_epid_tree(addOutputName);

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
        // ANCHOR FEMC cluster loop variables:
        // for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
        //     if(verbosity>1) std::cout << "\tFEMC:  cluster " << iclus << "\twith E = " << (_clusters_calo[kV1][kFEMC].at(iclus)).cluster_E << " GeV" << std::endl;
        // }

        // ANCHOR MC particle loop variables:
        // float * _mcpart_ID[imc]
        // float * _mcpart_ID_parent[imc]
        // float * _mcpart_PDG[imc]
        // float * _mcpart_E [imc]
        // float * _mcpart_px[imc]
        // float * _mcpart_py[imc]
        // float * _mcpart_pz[imc]

//	disreconstruction();

//	trackmatchingstudies();

	electronpid(i);

        clearClusterVectors();


    } // event loop end

	
//    trackmatchingstudiesSave();    

    write_epid_tree();
    std::cout << "all done :)" << std::endl;


}
