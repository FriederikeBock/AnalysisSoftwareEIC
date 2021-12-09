#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include "treeProcessing.h"
#include "event_utils.cxx"
#include "caloheader.h"
#include "clusterizer.cxx"

#include "disreconstruction.cxx"

#include "caloresolutionhistos.cxx"
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


void recoTreeProcessingSimple(
    Double_t maxNEvent          = -1,
    TString inFile              = "",
    TString inFileGeometry      = "geometry.root",
    TString addOutputName       = "",
    double electron_energy      = 10.,
    double hadron_energy        = 100.,
    Int_t verbosity             = 0
 ){

    bool do_reclus              = true;
    bool doCalibration          = true;
    
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0;
    bool brokenProjections = false;

    /*
    // make output directory
    TString dateForOutput = ReturnDateStr();
    outputDir = Form("treeProcessingSimple/%s",addOutputName.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // switch projection layer to most appropriate
    if (brokenProjections){
      _useAlternateForProjections = brokenProjections;
      std::cout << "calorimeter projections not available, using projections in TTL for matching" << std::endl;
    } else { 
      std::cout << "using calorimeter projections for matching" << std::endl;
    }
    */
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

    for (Int_t c = 0; c < 12; c++){
      cout << str_calorimeter[c] << "\t" << caloEnabled[c] << endl; 
    }

    Long64_t nEntriesTree                 = tt_event->GetEntries();
    std::cout << "Number of events in tree: " << nEntriesTree << std::endl;
    if(maxNEvent>0 && maxNEvent<nEntriesTree){
        nEntriesTree = maxNEvent;
        std::cout << "Will only analyze first " << maxNEvent << " events in the tree..." << std::endl;
    }
    _doClusterECalibration    = doCalibration;
    if(_doClusterECalibration){
        std::cout << "clusters will be energy-corrected and subsequently smeared to meet testbeam constant term!" << std::endl;
    }

    // DIS RECONSTRUCTION
    SetBeamEnergies(electron_energy,hadron_energy);
    SetCrossingAngle(25.e-3);
    SetDISFilename(addOutputName);
    InitializeDISRecon();

    float missedScatteredElectrons = 0;
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
        }
        if (verbosity > 1){
          for(int imc=0; imc<_nMCPart; imc++){
            std::cout << "MC part: \t" << imc << "\t E = " << _mcpart_E[imc] << "\tEta: " << _mcpart_Eta[imc]<< "\tPhi: " << _mcpart_Phi[imc] << "\t PDG: " << _mcpart_PDG[imc]  << "\t BCid: " << _mcpart_BCID[imc] << "\t HEPMC mother PDG: " << _hepmcp_PDG[_hepmcp_m1[_mcpart_BCID[imc]-1]-1] <<  std::endl;
          }
          
          for (int ihmc=0; ihmc < _nHepmcp; ihmc++){
            TVector3 hepvec(_hepmcp_px[ihmc],_hepmcp_py[ihmc],_hepmcp_pz[ihmc]);
            std::cout << "HEP MC part: \t" << ihmc << "\t E = " << _hepmcp_E[ihmc] << "\tEta: " << hepvec.Eta()<< "\tPhi: " << hepvec.Phi() << "\t PDG: " << _hepmcp_PDG[ihmc]  << "\t status: " << _hepmcp_status[ihmc] << "\t mother bcid: " << _hepmcp_m1[ihmc] << "\t";
            if (_hepmcp_m1[ihmc] > 0)
              std::cout <<  "PDG mother: " << _hepmcp_PDG[ _hepmcp_m1[ihmc]-1 ] << std::endl;
            else 
              std::cout <<  "source particle " << std::endl;
          }
          
          for (Int_t icalo = 0; icalo < maxcalo; icalo++){
            if (caloEnabled[icalo])
              std::cout << "towers " <<   str_calorimeter[icalo].Data() << "\t" << ReturnMaxTowerCalo(icalo) << std::endl;
          }
        }
        
        Int_t nTowers[maxcalo] = {_nTowers_FHCAL, _nTowers_FEMC, _nTowers_DRCALO, _nTowers_EEMC, _nTowers_CEMC,
                                  _nTowers_EHCAL, _nTowers_HCALIN, _nTowers_HCALOUT, _nTowers_LFHCAL, _nTowers_EEMCG,
                                  _nTowers_BECAL  };


        // Reset array entries for true ID
        for(Int_t itrk=0; itrk<_nTracks; itrk++){
          _track_trueID[itrk] = GetCorrectMCArrayEntry(_track_trueID[itrk]);
          if(verbosity>2){
            std::cout << "\tTrack: track " << itrk << " \t source: "<<_track_source[itrk] << "\twith true ID " << _track_trueID[itrk] << "\t px = " << _track_px[itrk] << " GeV ";
            if (_track_trueID[itrk] >= 0) std::cout << "\t pdg = " << _mcpart_PDG[(Int_t)_track_trueID[itrk]] << std::endl;
          }
        }
        
        // obtain labels for different track sources
        prepareMCMatchInfo();
        if (trackID_ScatteredElectron == -1 ){
         missedScatteredElectrons += 1; 
          if (verbosity > 0) 
            std::cout << "==================================== NO TRACK FOR THE SCATTERED ELECTRON WAS FOUND! ===============================================" << endl;
        }
        // run clusterizers normal calos
        for (int cal = 0; cal < maxcalo; cal++){
          if(do_reclus && nTowers[cal] > 0 && caloEnabled[cal]){
            for (int algo = 0; algo < _active_algo; algo++){
              if(verbosity>1) cout << "clusterizing " << str_clusterizer[algo].Data() <<  " for " << str_calorimeter[cal].Data() << endl;
              runclusterizer(algo, cal, seedE[cal], aggE[cal], primaryTrackSource);
            }
          }
        }
        if((do_reclus) && verbosity>1) cout << "done with clusterization!" << endl;

      // load correct clusterizer output also for centrally stored clusters, where appropriate
        for (Int_t icalo = 0; icalo < maxcalo; icalo++){
          for (Int_t ialgo = 0; ialgo < maxAlgo; ialgo++){
            if (loadClusterizerInput(ialgo, icalo) && verbosity>2) std::cout << str_calorimeter[icalo].Data() << "\t" << ialgo << endl;
          }
        }

        // set clusters matched to respective track for direct accessing
        SetClustersMatchedToTracks();

        // *****************************************************************************************************************
        // *****************************************************************************************************************
        // **********************   PUT your functions HERE ****************************************************************
        // *****************************************************************************************************************
        // *****************************************************************************************************************
       
      /*
        for(int itrk = 0; itrk<_nTracks; itrk++) {
          cout << "track ID = " << _track_ID[itrk] << "  track true ID = " << _track_trueID[itrk] << endl;
        }
        for(int imcp = 0; imcp<_nMCPart; imcp++) {
          cout << "mcpart ID = " << _mcpart_ID[imcp]  << "  mcpart ID parent = " << _mcpart_ID_parent[imcp] << endl;
        }
      */
        disreconstruction();
        // cleanup of variables
        clearClusterVectors();
        clearMCRecMatchVectors();

    } // event loop end
    
    std::cout << "fraction of events w/o scattered electron: " << float(missedScatteredElectrons/_nEventsTree) << std::endl;
    
    // *****************************************************************************************************************
    // *****************************************************************************************************************
    // **********************   PUT your save functions HERE ***********************************************************
    // *****************************************************************************************************************
    // *****************************************************************************************************************
  
    WriteDISTree(_nEventsTree);
    std::cout << "all done :)" << std::endl;
}
