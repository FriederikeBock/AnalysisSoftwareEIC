#include "treeProcessing.h"
#include "jet_finder.cxx"
#include "caloheader.h"
#include "clusterizer.cxx"
// #include "clusterizerNxN.cxx"
// #include "clusterizerV3.cxx"
// #include "clusterizerXN.cxx"
// #include "clusterizerNxNFEMC.cxx"
// #include "clusterizerV3FEMC.cxx"

#include "jetresolutionhistos.cxx"
#include "resolutionhistos.cxx"
#include "clusterstudies.cxx"
#include "trackingefficiency.cxx"
#include "hitstudies.cxx"
#include "trackmatchingstudies.cxx"

void treeProcessing(
    TString fileName            = "/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_pTHard5.root",
    TString addOutputName       = "",
    bool do_NxNclusterizer      = true,
    bool do_V3clusterizer       = true,
    bool do_XNclusterizer       = true,
    bool do_NxNclusterizerFEMC  = true,
    bool do_V3clusterizerFEMC   = true,
    bool do_XNclusterizerFEMC   = true,
    bool do_jetfinding          = false,
    // Double_t maxNEvent = 1e5,
    Double_t maxNEvent          = -1,
    Int_t verbosity             = 0
){
    // make output directory
    TString dateForOutput                       = ReturnDateStr();
    outputDir 						                  = Form("treeProcessing/%s_%s",addOutputName.Data(), dateForOutput.Data());
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
    _nEventsTree=0;
    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {
        // load current event
        tt_event->GetEntry(i);
        _nEventsTree++;

        // processing progress info
        if(i>0 && i%(nEntriesTree/(20)) ==0) cout << "//processed " << 100*(i)/nEntriesTree << "%"  << endl;
        // if(verbosity>1) cout << "event " << i << endl;

        // run clusterizers FHCAL
        if(do_NxNclusterizer){
            _do_NxNclusterizer = true;
            runclusterizer(kNxN, kFHCAL,0.5, 0.1,
                _nclusters_NxN_FHCAL,
                _clusters_NxN_FHCAL_E,
                _clusters_NxN_FHCAL_Eta,
                _clusters_NxN_FHCAL_Phi,
                _clusters_NxN_FHCAL_M02,
                _clusters_NxN_FHCAL_M20,
                _clusters_NxN_FHCAL_isMatched,
                _clusters_NxN_FHCAL_NTower,
                _clusters_NxN_FHCAL_trueID,
                _clusters_NxN_FHCAL_NtrueID,
                _clusters_NxN_FHCAL_X,
                _clusters_NxN_FHCAL_Y,
                _clusters_NxN_FHCAL_Z);
        }
        if(do_V3clusterizer){
            _do_V3clusterizer = true;
            runclusterizer(kV3, kFHCAL,0.5, 0.1,
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
        if(do_XNclusterizer){
            _do_XNclusterizer = true;
            runclusterizer(kXN, kFHCAL,0.5, 0.1,
                _nclusters_XN_FHCAL,
                _clusters_XN_FHCAL_E,
                _clusters_XN_FHCAL_Eta,
                _clusters_XN_FHCAL_Phi,
                _clusters_XN_FHCAL_M02,
                _clusters_XN_FHCAL_M20,
                _clusters_XN_FHCAL_isMatched,
                _clusters_XN_FHCAL_NTower,
                _clusters_XN_FHCAL_trueID,
                _clusters_XN_FHCAL_NtrueID,
                _clusters_XN_FHCAL_X,
                _clusters_XN_FHCAL_Y,
                _clusters_XN_FHCAL_Z);
        }
        // run clusterizers FEMC
        if(do_NxNclusterizerFEMC){
            _do_NxNclusterizerFEMC = true;
            runclusterizer(kNxN, kFEMC,0.5, 0.1,
                _nclusters_NxN_FEMC,
                _clusters_NxN_FEMC_E,
                _clusters_NxN_FEMC_Eta,
                _clusters_NxN_FEMC_Phi,
                _clusters_NxN_FEMC_M02,
                _clusters_NxN_FEMC_M20,
                _clusters_NxN_FEMC_isMatched,
                _clusters_NxN_FEMC_NTower,
                _clusters_NxN_FEMC_trueID,
                _clusters_NxN_FEMC_NtrueID,
                _clusters_NxN_FEMC_X,
                _clusters_NxN_FEMC_Y,
                _clusters_NxN_FEMC_Z);
        }
        if(do_V3clusterizerFEMC){
            _do_V3clusterizerFEMC = true;
            runclusterizer(kV3, kFEMC,0.5, 0.1,
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
        if(do_XNclusterizerFEMC){
            _do_XNclusterizerFEMC = true;
            runclusterizer(kXN, kFEMC,0.5, 0.1,
                _nclusters_XN_FEMC,
                _clusters_XN_FEMC_E,
                _clusters_XN_FEMC_Eta,
                _clusters_XN_FEMC_Phi,
                _clusters_XN_FEMC_M02,
                _clusters_XN_FEMC_M20,
                _clusters_XN_FEMC_isMatched,
                _clusters_XN_FEMC_NTower,
                _clusters_XN_FEMC_trueID,
                _clusters_XN_FEMC_NtrueID,
                _clusters_XN_FEMC_X,
                _clusters_XN_FEMC_Y,
                _clusters_XN_FEMC_Z);
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
        for(Int_t itrk=0; itrk<_nTracks; itrk++){
            if(verbosity>1) cout << "\tTrack: track " << itrk << "\twith true ID " << _track_trueID[itrk]-1 << "\tand X = " << _track_px[itrk] << " cm" << endl;

            if(_do_jetfinding){
                TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
                if(trackvec.Eta()<0) continue;
                float Etrack = TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)+TMath::Power(_track_pz[itrk],2));
                // create track vector for jet finder
                jetf_track_px.push_back(_track_px[itrk]);
                jetf_track_py.push_back(_track_py[itrk]);
                jetf_track_pz.push_back(_track_pz[itrk]);
                jetf_track_E.push_back(Etrack);
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
        if(do_V3clusterizer){
            for(Int_t iclus=0; iclus<_nclusters_V3_FHCAL; iclus++){
                if(!_clusters_V3_FHCAL_isMatched[iclus]){
                    double pt = _clusters_V3_FHCAL_E[iclus] / cosh(_clusters_V3_FHCAL_Eta[iclus]);
                    double px = pt * cos(_clusters_V3_FHCAL_Phi[iclus]);
                    double py = pt * sin(_clusters_V3_FHCAL_Phi[iclus]);
                    double pz = pt * sinh(_clusters_V3_FHCAL_Eta[iclus]);
                    jetf_all_px.push_back(px);
                    jetf_all_py.push_back(py);
                    jetf_all_pz.push_back(pz);
                    jetf_all_E.push_back(_clusters_V3_FHCAL_E[iclus]);
                    jetf_hcal_px.push_back(px);
                    jetf_hcal_py.push_back(py);
                    jetf_hcal_pz.push_back(pz);
                    jetf_hcal_E.push_back(_clusters_V3_FHCAL_E[iclus]);
                    jetf_calo_px.push_back(px);
                    jetf_calo_py.push_back(py);
                    jetf_calo_pz.push_back(pz);
                    jetf_calo_E.push_back(_clusters_V3_FHCAL_E[iclus]);
                    jetf_all_px.push_back(px);
                    jetf_all_py.push_back(py);
                    jetf_all_pz.push_back(pz);
                    jetf_all_E.push_back(_clusters_V3_FHCAL_E[iclus]);
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
        if(do_V3clusterizerFEMC){
            for(Int_t iclus=0; iclus<_nclusters_V3_FEMC; iclus++){
                if(!_clusters_V3_FEMC_isMatched[iclus]){
                    double pt = _clusters_V3_FEMC_E[iclus] / cosh(_clusters_V3_FEMC_Eta[iclus]);
                    double px = pt * cos(_clusters_V3_FEMC_Phi[iclus]);
                    double py = pt * sin(_clusters_V3_FEMC_Phi[iclus]);
                    double pz = pt * sinh(_clusters_V3_FEMC_Eta[iclus]);
                    jetf_full_px.push_back(px);
                    jetf_full_py.push_back(py);
                    jetf_full_pz.push_back(pz);
                    jetf_full_E.push_back(_clusters_V3_FEMC_E[iclus]);
                    jetf_calo_px.push_back(px);
                    jetf_calo_py.push_back(py);
                    jetf_calo_pz.push_back(pz);
                    jetf_calo_E.push_back(_clusters_V3_FEMC_E[iclus]);
                    jetf_all_px.push_back(px);
                    jetf_all_py.push_back(py);
                    jetf_all_pz.push_back(pz);
                    jetf_all_E.push_back(_clusters_V3_FEMC_E[iclus]);
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
                if(truevec.Eta()<0) continue;
                if(verbosity>1) cout << "\tMC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;
                //  create truth vector for jet finder
                jetf_truth_px.push_back(_mcpart_px[imc]);
                jetf_truth_py.push_back(_mcpart_py[imc]);
                jetf_truth_pz.push_back(_mcpart_pz[imc]);
                jetf_truth_E.push_back(_mcpart_E[imc]);
                // fill charged only inputs (reject neutral particles)
                if(abs(_mcpart_PDG[imc])==22) continue; // photon
                if(abs(_mcpart_PDG[imc])==2112) continue; // neutron
                jetf_truthcharged_px.push_back(_mcpart_px[imc]);
                jetf_truthcharged_py.push_back(_mcpart_py[imc]);
                jetf_truthcharged_pz.push_back(_mcpart_pz[imc]);
                jetf_truthcharged_E.push_back(_mcpart_E[imc]);
            }

            // ANCHOR JET FINDING
            // truth jets
            auto jetsTrue = findJets(1.0, "anti-kt", jetf_truth_px, jetf_truth_py, jetf_truth_pz, jetf_truth_E);
            if(verbosity>1)cout << "found " << std::get<1>(jetsTrue).size() << " true jets" << endl;        // printJets(std::get<1>(jetsTrue));

            auto jetsTrueCharged = findJets(1.0, "anti-kt", jetf_truthcharged_px, jetf_truthcharged_py, jetf_truthcharged_pz, jetf_truthcharged_E);
            if(verbosity>1)cout << "found " << std::get<1>(jetsTrueCharged).size() << " true charged jets" << endl;        // printJets(std::get<1>(jetsTrue));

            // track-based jets (rec)
            auto jetsTrackRec = findJets(1.0, "anti-kt", jetf_track_px, jetf_track_py, jetf_track_pz, jetf_track_E);
            if(verbosity>1)cout << "found " << std::get<1>(jetsTrackRec).size() << " rec track jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // full jets (rec)
            auto jetsFullRec = findJets(1.0, "anti-kt", jetf_full_px, jetf_full_py, jetf_full_pz, jetf_full_E);
            if(verbosity>1) cout << "found " << std::get<1>(jetsFullRec).size() << " rec full jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // hcal jets (rec)
            auto jetsHcalRec = findJets(1.0, "anti-kt", jetf_hcal_px, jetf_hcal_py, jetf_hcal_pz, jetf_hcal_E);
            if(verbosity>1) cout << "found " << std::get<1>(jetsHcalRec).size() << " rec hcal jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // calo jets (rec)
            auto jetsCaloRec = findJets(1.0, "anti-kt", jetf_calo_px, jetf_calo_py, jetf_calo_pz, jetf_calo_E);
            if(verbosity>1) cout << "found " << std::get<1>(jetsCaloRec).size() << " rec calo jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            // all jets (rec)
            auto jetsAllRec = findJets(1.0, "anti-kt", jetf_all_px, jetf_all_py, jetf_all_pz, jetf_all_E);
            if(verbosity>1) cout << "found " << std::get<1>(jetsAllRec).size() << " rec all jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

            jetresolutionhistos(jetsTrackRec,jetsTrueCharged,0);
            jetresolutionhistos(jetsFullRec,jetsTrue,1);
            jetresolutionhistos(jetsHcalRec,jetsTrue,2);
            jetresolutionhistos(jetsCaloRec,jetsTrue,3);
            jetresolutionhistos(jetsAllRec,jetsTrue,4);
// TString jettype[njettypes] = {"track", "full","hcal","calo","all"};

        }
        resolutionhistos();
        clusterstudies();
        // if(_do_NxNclusterizer) trackmatchingstudies(kNxN, kFHCAL,true);
        // if(_do_NxNclusterizerFEMC) trackmatchingstudies(kNxN, kFEMC,true);
        // if(_do_V3clusterizer) trackmatchingstudies(kV3, kFHCAL,true);
        // if(_do_V3clusterizerFEMC) trackmatchingstudies(kV3, kFEMC,true);

    } // event loop end
    jetresolutionhistosSave();
    resolutionhistosSave();
    clusterstudiesSave();
    trackingefficiencyhistosSave();
    hitstudiesSave();
    // trackmatchingstudiesSave();
    clusterizerSave();
}
