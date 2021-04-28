#include "../common/binningheader.h"
#include "treeProcessing.h"
#include "jet_finder.cxx"
#include "caloheader.h"
#include "clusterizer.cxx"

#include "jetresolutionhistos.cxx"
#include "resolutionhistos.cxx"
#include "clusterstudies.cxx"
#include "trackingefficiency.cxx"
#include "hitstudies.cxx"
#include "trackmatchingstudies.cxx"

void treeProcessing(
    TString fileName            = "/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_pTHard5.root",
    TString addOutputName       = "",
    bool do_3x3clusterizer      = true,
    bool do_5x5clusterizer      = true,
    bool do_V3clusterizer       = true,
    bool do_MAclusterizer       = true,
    bool do_C3clusterizer       = true,
    bool do_C5clusterizer       = true,
    bool do_3x3clusterizerFEMC  = true,
    bool do_5x5clusterizerFEMC  = true,
    bool do_V3clusterizerFEMC   = true,
    bool do_MAclusterizerFEMC   = true,
    bool do_C3clusterizerFEMC   = true,
    bool do_C5clusterizerFEMC   = true,
    bool do_jetfinding          = false,
    // Double_t maxNEvent = 1e5,
    bool hasTiming              = true,
    bool isALLSILICON           = true, 
    Double_t maxNEvent          = -1,
    Int_t verbosity             = 0,
    int granularity_factor      = 1,
    bool doCalibration          = false
    // Defaults to tracking from all layers.
    unsigned int primaryTrackSource = 0
){
    // make output directory
    TString dateForOutput                       = ReturnDateStr();
    outputDir 						                  = Form("treeProcessing/%s",addOutputName.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    if(do_jetfinding) _do_jetfinding = true;

    // load tree
    TTree *const tt_event =  (TTree *) (new TFile(fileName.Data(), "READ"))->Get("event_tree");
    if(!tt_event){ cout << "tree not found... returning!"<< endl; return;}

    // load all branches (see header)
    SetBranchAddressesTree(tt_event);

    Long64_t nEntriesTree                 = tt_event->GetEntries();
    cout << "Number of events in tree: " << nEntriesTree << endl;
    if(maxNEvent>0 && maxNEvent<nEntriesTree){
        nEntriesTree = maxNEvent;
        cout << "Will only analyze first " << maxNEvent << " events in the tree..." << endl;
    }
    _do_TimingStudies         = hasTiming;
    _is_ALLSILICON            = isALLSILICON;
    _granularityCalo          = granularity_factor;
    _doClusterECalibration    = doCalibration;
    if(_doClusterECalibration){
        cout << "clusters will be energy-corrected and subsequently smeared to meet testbeam constant term!" << endl;
    }
    _nEventsTree=0;
    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {
        // load current event
        tt_event->GetEntry(i);
        _nEventsTree++;

        // processing progress info
        if(i>0 && nEntriesTree>100 && i%(nEntriesTree/(20))==0) cout << "//processed " << 100*(i)/nEntriesTree << "%"  << endl;
        if(verbosity>1) cout << "event " << i << endl;

        // calculate useful quantities
        for(Int_t imc=0; imc<_nMCPart; imc++){
            TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
            _mcpart_Eta[imc]=truevec.Eta();
        }
        float seed_E = 0.5;
        float aggregation_E = 0.1;
        // run clusterizers FHCAL
        if(do_3x3clusterizer){
            _do_3x3clusterizer = true;
            if(verbosity>1) cout << "clusterizing 3x3 for FHCAL" << endl;
            runclusterizer(k3x3, kFHCAL,seed_E, aggregation_E,
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
        if(do_5x5clusterizer){
            _do_5x5clusterizer = true;
            if(verbosity>1) cout << "clusterizing 5x5 for FHCAL" << endl;
            runclusterizer(k5x5, kFHCAL,seed_E, aggregation_E,
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
        if(do_V3clusterizer){
            _do_V3clusterizer = true;
            if(verbosity>1) cout << "clusterizing V3 for FHCAL" << endl;
            runclusterizer(kV3, kFHCAL,seed_E, aggregation_E,
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
        if(do_MAclusterizer){
            _do_MAclusterizer = true;
            if(verbosity>1) cout << "clusterizing MA for FHCAL" << endl;
            runclusterizer(kMA, kFHCAL,seed_E, aggregation_E,
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
        if(do_C3clusterizer){
            _do_C3clusterizer = true;
            if(verbosity>1) cout << "clusterizing C3 for FHCAL" << endl;
            runclusterizer(kC3, kFHCAL,seed_E, aggregation_E,
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
        if(do_C5clusterizer){
            _do_C5clusterizer = true;
            if(verbosity>1) cout << "clusterizing C5 for FHCAL" << endl;
            runclusterizer(kC5, kFHCAL,seed_E, aggregation_E,
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
        // run clusterizers FEMC
        if(do_3x3clusterizerFEMC){
            _do_3x3clusterizerFEMC = true;
            if(verbosity>1) cout << "clusterizing 3x3 for FEMC" << endl;
            runclusterizer(k3x3, kFEMC,seed_E, aggregation_E,
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
        if(do_5x5clusterizerFEMC){
            _do_5x5clusterizerFEMC = true;
            if(verbosity>1) cout << "clusterizing 5x5 for FEMC" << endl;
            runclusterizer(k5x5, kFEMC,seed_E, aggregation_E,
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
        if(do_V3clusterizerFEMC){
            _do_V3clusterizerFEMC = true;
            if(verbosity>1) cout << "clusterizing V3 for FEMC" << endl;
            runclusterizer(kV3, kFEMC,seed_E, aggregation_E,
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
        if(do_MAclusterizerFEMC){
            _do_MAclusterizerFEMC = true;
            if(verbosity>1) cout << "clusterizing MA for FEMC" << endl;
            runclusterizer(kMA, kFEMC,seed_E, aggregation_E,
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
        if(do_C3clusterizerFEMC){
            _do_C3clusterizerFEMC = true;
            if(verbosity>1) cout << "clusterizing C3 for FEMC" << endl;
            runclusterizer(kC3, kFEMC,seed_E, aggregation_E,
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
        if(do_C5clusterizerFEMC){
            _do_C5clusterizerFEMC = true;
            if(verbosity>1) cout << "clusterizing C5 for FEMC" << endl;
            runclusterizer(kC5, kFEMC,seed_E, aggregation_E,
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
        if(0){ //do_V1clusterizerDRCALO
            // _do_V1clusterizerDRCALO = true;
            float seed_E_DRCALO = 0.2;
            float aggregation_E_DRCALO = 0.001;
            if(verbosity>1) cout << "clusterizing V1 for DRCALO" << endl;
            runclusterizer(kV1, kDRCALO,seed_E_DRCALO, aggregation_E_DRCALO,
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
        }

        // ANCHOR Hits loop variables:
        // float* _hits_x[ihit]
        // float* _hits_y[ihit]
        // float* _hits_z[ihit]
        // float* _hits_t[ihit]
        // for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
        //     if(verbosity>1) cout << "\tHIT: hit " << ihit << "\tin layer " << _hits_layerID[ihit] << "\twith X = " << _hits_x[ihit] << " cm" << endl;
        // }
        hitstudies();

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
            if(verbosity>1) cout << "\tTrack: track " << itrk << "\twith true ID " << _track_trueID[itrk]-1 << "\tand X = " << _track_px[itrk] << " cm" << endl;

            // Skip tracks that aren't selected.
            if(_track_source[itrk] != primaryTrackSource) {
                continue;
            }

            if(_do_jetfinding){
                TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
                float Etrack = TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)+TMath::Power(_track_pz[itrk],2) + TMath::Power(massHypothesis, 2));
                // create track vector for jet finder
                jetf_track_px.push_back(_track_px[itrk]);
                jetf_track_py.push_back(_track_py[itrk]);
                jetf_track_pz.push_back(_track_pz[itrk]);
                jetf_track_E.push_back(Etrack);
                if(trackvec.Eta()<0) continue;
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
        trackingefficiency();
        trackingresolution();
        trackingcomparison();

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
        //     if(verbosity>1) cout << "\tProjection: proj " << iproj << "\tin layer " << _track_ProjLayer[iproj] << "\twith X = " << _track_Proj_x[iproj] << " cm" << endl;
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
        std::vector<float> jetf_calo_E;
        std::vector<float> jetf_calo_px;
        std::vector<float> jetf_calo_py;
        std::vector<float> jetf_calo_pz;
        // for(Int_t iclus=0; iclus<_nclusters_FHCAL; iclus++){
        //     if(verbosity>1) cout << "\tcls " << iclus << "\tE " << _clusters_FHCAL_E[iclus] << "\tEta " << _clusters_FHCAL_Eta[iclus] << "\tPhi " << _clusters_FHCAL_Phi[iclus] << "\tntowers: " << _clusters_FHCAL_NTower[iclus] << "\ttrueID: " << _clusters_FHCAL_trueID[iclus] << endl;
        // }
                // apply calibration if desired
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
        if(do_MAclusterizer){
            for(Int_t iclus=0; iclus<_nclusters_MA_FHCAL; iclus++){
                if(!_clusters_MA_FHCAL_isMatched[iclus]){
                    double pt = _clusters_MA_FHCAL_E[iclus] / cosh(_clusters_MA_FHCAL_Eta[iclus]);
                    double px = pt * cos(_clusters_MA_FHCAL_Phi[iclus]);
                    double py = pt * sin(_clusters_MA_FHCAL_Phi[iclus]);
                    double pz = pt * sinh(_clusters_MA_FHCAL_Eta[iclus]);
                    jetf_all_px.push_back(px);
                    jetf_all_py.push_back(py);
                    jetf_all_pz.push_back(pz);
                    jetf_all_E.push_back(_clusters_MA_FHCAL_E[iclus]);
                    jetf_hcal_px.push_back(px);
                    jetf_hcal_py.push_back(py);
                    jetf_hcal_pz.push_back(pz);
                    jetf_hcal_E.push_back(_clusters_MA_FHCAL_E[iclus]);
                    jetf_calo_px.push_back(px);
                    jetf_calo_py.push_back(py);
                    jetf_calo_pz.push_back(pz);
                    jetf_calo_E.push_back(_clusters_MA_FHCAL_E[iclus]);
                    jetf_all_px.push_back(px);
                    jetf_all_py.push_back(py);
                    jetf_all_pz.push_back(pz);
                    jetf_all_E.push_back(_clusters_MA_FHCAL_E[iclus]);
                }
            }
        }

        // ANCHOR FEMC cluster loop variables:
        // float* _clusters_FEMC_E[iclus];
        // float* _clusters_FEMC_Eta[iclus];
        // float* _clusters_FEMC_Phi[iclus];
        // int* _clusters_FEMC_NTower[iclus];
        // float* _clusters_FEMC_trueID[iclus];
        // for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
        //     if(verbosity>1) cout << "\tFEMC:  cluster " << iclus << "\twith E = " << _clusters_FEMC_E[iclus] << " GeV" << endl;
        // }
        if(do_MAclusterizerFEMC){
            for(Int_t iclus=0; iclus<_nclusters_MA_FEMC; iclus++){
                if(!_clusters_MA_FEMC_isMatched[iclus]){
                    double pt = _clusters_MA_FEMC_E[iclus] / cosh(_clusters_MA_FEMC_Eta[iclus]);
                    double px = pt * cos(_clusters_MA_FEMC_Phi[iclus]);
                    double py = pt * sin(_clusters_MA_FEMC_Phi[iclus]);
                    double pz = pt * sinh(_clusters_MA_FEMC_Eta[iclus]);
                    jetf_full_px.push_back(px);
                    jetf_full_py.push_back(py);
                    jetf_full_pz.push_back(pz);
                    jetf_full_E.push_back(_clusters_MA_FEMC_E[iclus]);
                    jetf_calo_px.push_back(px);
                    jetf_calo_py.push_back(py);
                    jetf_calo_pz.push_back(pz);
                    jetf_calo_E.push_back(_clusters_MA_FEMC_E[iclus]);
                    jetf_all_px.push_back(px);
                    jetf_all_py.push_back(py);
                    jetf_all_pz.push_back(pz);
                    jetf_all_E.push_back(_clusters_MA_FEMC_E[iclus]);
                }
            }
        }
        // ANCHOR FHCAL tower loop variables:
        // float* _tower_FHCAL_E[itwrH]
        // int* _tower_FHCAL_iEta[itwrH]
        // int* _tower_FHCAL_iPhi[itwrH]
        // float* _tower_FHCAL_trueID[itwrH]
        // for(Int_t itwrH=0; itwrH<_nTowers_FHCAL; itwrH++){
        //     if(verbosity>1) cout << "\tFHCAL: tower " << itwrH << "\twith E = " << _tower_FHCAL_E[itwrH] << " GeV" << endl;
        // }

        // ANCHOR FEMC cluster loop variables:
        // float* _tower_FEMC_E[itwrH]
        // int* _tower_FEMC_iEta[itwrH]
        // int* _tower_FEMC_iPhi[itwrH]
        // float* _tower_FEMC_trueID[itwrH]
        // for(Int_t itwrE=0; itwrE<_nTowers_FEMC; itwrE++){
        //     if(verbosity>1) cout << "\tFEMC:  tower " << itwrE << "\twith E = " << _tower_FEMC_E[itwrE] << " GeV" << endl;
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
                    if(verbosity>1) cout << "\tMC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;
                    //  create truth vector for jet finder
                    jetf_truth_px.push_back(_mcpart_px[imc]);
                    jetf_truth_py.push_back(_mcpart_py[imc]);
                    jetf_truth_pz.push_back(_mcpart_pz[imc]);
                    jetf_truth_E.push_back(_mcpart_E[imc]);
//                 } 
                // fill charged only inputs (reject neutral particles)
//                 cout << _mcpart_PDG[imc] << endl;
                if(abs(_mcpart_PDG[imc])==211 || abs(_mcpart_PDG[imc])==321 || abs(_mcpart_PDG[imc])==2212 || abs(_mcpart_PDG[imc])==11 || abs(_mcpart_PDG[imc])==13){ 
                    jetf_truthcharged_px.push_back(_mcpart_px[imc]);
                    jetf_truthcharged_py.push_back(_mcpart_py[imc]);
                    jetf_truthcharged_pz.push_back(_mcpart_pz[imc]);
                    jetf_truthcharged_E.push_back(_mcpart_E[imc]);
                }
            }

            // ANCHOR JET FINDING
            // truth jets
            auto jetsTrue = findJets(0.5, "anti-kt", jetf_truth_px, jetf_truth_py, jetf_truth_pz, jetf_truth_E);
            if(verbosity>1)cout << "found " << std::get<1>(jetsTrue).size() << " true jets" << endl;        // printJets(std::get<1>(jetsTrue));

            auto jetsTrueCharged = findJets(0.5, "anti-kt", jetf_truthcharged_px, jetf_truthcharged_py, jetf_truthcharged_pz, jetf_truthcharged_E);
            if(verbosity>1)cout << "found " << std::get<1>(jetsTrueCharged).size() << " true charged jets" << endl;        // printJets(std::get<1>(jetsTrue));

            // track-based jets (rec)
            auto jetsTrackRec = findJets(0.5, "anti-kt", jetf_track_px, jetf_track_py, jetf_track_pz, jetf_track_E);
            if(verbosity>1)cout << "found " << std::get<1>(jetsTrackRec).size() << " rec track jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // full jets (rec)
            auto jetsFullRec = findJets(0.5, "anti-kt", jetf_full_px, jetf_full_py, jetf_full_pz, jetf_full_E);
            if(verbosity>1) cout << "found " << std::get<1>(jetsFullRec).size() << " rec full jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // hcal jets (rec)
            // auto jetsHcalRec = findJets(0.5, "anti-kt", jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz, jetf_hcal_E);
            // if(verbosity>1) cout << "found " << std::get<1>(jetsHcalRec).size() << " rec hcal jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // calo jets (rec)
            auto jetsCaloRec = findJets(0.5, "anti-kt", jetf_calo_px, jetf_calo_py, jetf_calo_pz, jetf_calo_E);
            if(verbosity>1) cout << "found " << std::get<1>(jetsCaloRec).size() << " rec calo jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // all jets (rec)
            auto jetsAllRec = findJets(0.5, "anti-kt", jetf_all_px, jetf_all_py, jetf_all_pz, jetf_all_E);
            if(verbosity>1) cout << "found " << std::get<1>(jetsAllRec).size() << " rec all jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            jetresolutionhistos(jetsTrackRec,jetsTrueCharged,0);
            jetresolutionhistos(jetsFullRec,jetsTrue,1);
            // jetresolutionhistos(jetsHcalRec,jetsTrue,2);
            jetresolutionhistos(jetsCaloRec,jetsTrue,3);
            jetresolutionhistos(jetsAllRec,jetsTrue,4);
// TString jettype[njettypes] = {"track", "full","hcal","calo","all"};
        }
        if(verbosity>1) cout << "running clusterstudies" << endl;
        clusterstudies();
        if(verbosity>1) cout << "running resolutionhistos" << endl;
        resolutionhistos();
        if(verbosity>1) cout << "loop done ... next event" << endl;
        trackmatchingstudies();

    } // event loop end
    cout << "running jetresolutionhistosSave" << endl;
    jetresolutionhistosSave();
    cout << "running resolutionhistosSave" << endl;
    resolutionhistosSave();
    cout << "running clusterstudiesSave" << endl;
    clusterstudiesSave();
    cout << "running trackingefficiencyhistosSave" << endl;
    trackingefficiencyhistosSave();
    cout << "running trackingresolutionhistosSave" << endl;
    trackingresolutionhistosSave();
    cout << "running trackingcomparisonhistosSave" << endl;
    trackingcomparisonhistosSave();
    cout << "running hitstudiesSave" << endl;
    hitstudiesSave();
    cout << "running trackmatchingstudiesSave" << endl;
    trackmatchingstudiesSave();
    cout << "running clusterizerSave" << endl;
    clusterizerSave();
    cout << "all done :)" << endl;
}
