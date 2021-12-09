#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMath.h>
#include <TVector3.h>
#include <iostream>
#include <fstream>

#include "disreconstruction.cxx"

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

void recoTreeProcessing(
    Double_t maxNEvent		= -1,
    TString inFile              = "",
    TString inFileGeometry      = "geometry.root",
    TString addOutputName       = "",
    double electron_energy	= 10.,
    double hadron_energy 	= 100.
 ){

    bool do_reclus              = true;
    bool do_jetfinding          = false;
    bool runCaloRes             = true,
    bool doCalibration          = true,
    Int_t verbosity             = 0;
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0;
    std::string jetAlgorithm    = "anti-kt";
    double jetR                 = 0.5;
    double tracked_jet_max_pT   = 30;
    bool removeTracklets        = false

/*
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
*/
	
    TChain *const tt_event = new TChain("event_tree");
    tt_event->Add(inFile);



    // // load geometry tree
    tt_geometry =  (TTree *) (new TFile(inFileGeometry.Data(), "READ"))->Get("geometry_tree");
    if(!tt_geometry){ std::cout << "geometry tree not found... returning!"<< std::endl; return;}
    // load all branches (see header)
    SetBranchAddressesTree(tt_event);
    SetBranchAddressesGeometryTree(tt_geometry);
    SetGeometryIndices();

    for (Int_t c = 0; c < 12; c++){
      std::cout << str_calorimeter[c] << "\t" << caloEnabled[c] << std::endl; 
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

/*
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
*/

	// DIS RECONSTRUCTION
	SetBeamEnergies(electron_energy,hadron_energy);
	SetCrossingAngle(25.e-3);
	SetDISFilename(addOutputName);
	InitializeDISRecon();

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
              if(verbosity>1) std::cout << "clusterizing " << str_clusterizer[algo].Data() <<  " for " << str_calorimeter[cal].Data() << std::endl;
              runclusterizer(algo, cal, seedE[cal], aggE[cal], primaryTrackSource);
            }
          }
        }

        if(do_reclus && _nTowers_FOCAL && caloEnabled[kFOCAL]){
            float seed_E_FOCAL = 0.1;
            float aggregation_E_FOCAL = 0.001;
            if(kMA<_active_algo){
                if(verbosity>1) std::cout << "clusterizing MA for FOCAL" << std::endl;
                runclusterizer(kMA, kFOCAL,seed_E_FOCAL, aggregation_E_FOCAL, primaryTrackSource);
            }
        } 

        if(do_reclus && _nTowers_DRCALO && caloEnabled[kDRCALO]){ //do_V1clusterizerDRCALO
          if(verbosity>1) std::cout << "clusterizing V1 for DRCALO" << std::endl;
          runclusterizer(kV1, kDRCALO, seedE[kDRCALO], aggE[kDRCALO], primaryTrackSource);
        }
        if((do_reclus) && verbosity>1) std::cout << "done with clusterization!" << std::endl;

        double massHypothesis = 0.139;
        for(Int_t itrk=0; itrk<_nTracks; itrk++){
            _track_trueID[itrk] = GetCorrectMCArrayEntry(_track_trueID[itrk]);
            // if(verbosity>1) std::cout << "\tTrack: track " << itrk << "\twith true ID " << _track_trueID[itrk] << "\tand X = " << _track_px[itrk] << " cm" << std::endl;
        }
        // obtain labels for different track sources
        prepareMCMatchInfo();

        int fwdCaloID = kFHCAL;
        if (!caloEnabled[kFHCAL]) fwdCaloID = kLFHCAL;

        for (Int_t icalo = 0; icalo < maxcalo; icalo++){
          for (Int_t ialgo = 0; ialgo < maxAlgo; ialgo++){
            if (loadClusterizerInput(ialgo, icalo) && verbosity>2) std::cout << str_calorimeter[icalo].Data() << "\t" << ialgo << std::endl;
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

	disreconstruction();
        clearClusterVectors();
        clearMCRecMatchVectors();

    } // event loop end

    WriteDISTree();

}



// MAIN FUNCTION for non-ROOT compilation (prototype)
int main( int argc, char* argv[] )
{

    // Default arguments for treeProcessing
    TString inFile              = "";
    TString inFileGeometry      = "geometry.root";
    TString addOutputName       = "";
    Double_t maxNEvent          = -1;
    bool do_reclus              = true;
    bool do_jetfinding          = false;
    bool runCaloRes             = true;
    bool doCalibration          = false;
    Int_t verbosity             = 0;
    // Defaults to tracking from all layers.
    unsigned short primaryTrackSource = 0;
    std::string jetAlgorithm    = "anti-kt";
    double jetR                 = 0.5;
    double tracked_jet_max_pT   = 30;
    bool removeTracklets        = false;

    // Import main call arguments
        TString import;
        if( argc >  1 ) inFile                 = argv[1];
        if( argc >  2 ) inFileGeometry                  = argv[2];
        if( argc >  3 ) addOutputName          = argv[3];
        if( argc >  4 ) { // numberOfBins
            std::istringstream sstr(argv[4]);
            sstr >> maxNEvent;
        }
        if( argc >  5 ) { // UseTHnSparse
            import = argv[5];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  do_reclus = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  do_reclus = kFALSE;
        }
        if( argc >  6 ) { // UseTHnSparse
            import = argv[6];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  do_jetfinding = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  do_jetfinding = kFALSE;
        }
        if( argc >  7 ) { // UseTHnSparse
            import = argv[7];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  runCaloRes = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  runCaloRes = kFALSE;
        }
        if( argc >  8 ) { // UseTHnSparse
            import = argv[8];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  doCalibration = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  doCalibration = kFALSE;
        }
        if( argc >  9 ) {
          std::istringstream sstr(argv[9]);
          sstr >> verbosity;
        }
        if( argc > 10 ) {
          std::istringstream sstr(argv[10]);
          sstr >> primaryTrackSource;
        }
        if( argc > 11 ) jetAlgorithm = argv[11];
        if( argc > 12 ) { // numberOfBins
            std::stringstream sstr(argv[12]);
            sstr >> jetR;
        } if( argc > 13 ) { // numberOfBins
            std::stringstream sstr(argv[13]);
            sstr >> tracked_jet_max_pT;
        } if( argc > 14 ) { // UseTHnSparse
            import = argv[14];
            if( import.EqualTo("kTRUE") || import.EqualTo("true") )  removeTracklets = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("false") )  removeTracklets = kFALSE;
        }

    // Function call ExtractSignalV2
        treeProcessing(
            inFile,
            inFileGeometry,
            addOutputName,
            maxNEvent,
            do_reclus,
            do_jetfinding,
            runCaloRes,
            doCalibration,
            verbosity,
            primaryTrackSource,
            jetAlgorithm,
            jetR,
            tracked_jet_max_pT,
            removeTracklets
        );
    return 0;


}
