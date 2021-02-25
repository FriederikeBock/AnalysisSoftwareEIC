#include "jet_finder.cxx"
#include "treeProcessing.h"
#include "resolutionhistos.cxx"


void treeProcessing(
    TString fileName                = "/media/nschmidt/local/EIC_running/Modular_ALLSILICON-FTTLS3LC-ETTL-CTTL_testing/eventtree.root",
    Double_t maxNEvent = -1,
    Int_t verbosity = 1
){
    // make output directory
    gSystem->Exec("mkdir -p treeProcessing");

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

    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {
        // load current event
        tt_event->GetEntry(i);

        // processing progress info
        if(i>0 && i%(nEntriesTree/(10)) ==0) cout << "//processed " << 100*(i)/nEntriesTree << "%"  << endl;

        // ANCHOR Hits loop variables:
        // float* _hits_x[ihit]
        // float* _hits_y[ihit]
        // float* _hits_z[ihit]
        // float* _hits_t[ihit]
        for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
            if(verbosity>1) cout << "\tHIT: hit " << ihit << "\tin layer " << _hits_layerID[ihit] << "\twith X = " << _hits_x[ihit] << " cm" << endl;
        }

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
        for(Int_t itrk=0; itrk<_nTracks; itrk++){
            if(verbosity>1) cout << "\tTrack: track " << itrk << "\twith true ID " << _track_trueID[itrk] << "\tand X = " << _track_px[itrk] << " cm" << endl;

            // create track vector for jet finder
            jetf_track_px.push_back(_track_px[itrk]);
            jetf_track_py.push_back(_track_py[itrk]);
            jetf_track_pz.push_back(_track_pz[itrk]);
            jetf_track_E.push_back(TMath::Sqrt(TMath::Power(_track_px[itrk],2)+TMath::Power(_track_py[itrk],2)+TMath::Power(_track_pz[itrk],2)));
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
        for(Int_t iproj=0; iproj<_nProjections; iproj++){
            if(verbosity>1) cout << "\tProjection: proj " << iproj << "\tin layer " << _track_ProjLayer[iproj] << "\twith X = " << _track_Proj_x[iproj] << " cm" << endl;
        }

        // ANCHOR FHCAL cluster loop variables:
        // float* _cluster_FHCAL_E[iclus];
        // float* _cluster_FHCAL_Eta[iclus];
        // float* _cluster_FHCAL_Phi[iclus];
        // int* _cluster_FHCAL_NTower[iclus];
        // float* _cluster_FHCAL_trueID[iclus];
        for(Int_t iclus=0; iclus<_nclusters_FHCAL; iclus++){
            if(verbosity>1) cout << "\tFHCAL: cluster " << iclus << "\twith E = " << _cluster_FHCAL_E[iclus] << " GeV" << endl;
        }

        // ANCHOR FEMC cluster loop variables:
        // float* _cluster_FEMC_E[iclus];
        // float* _cluster_FEMC_Eta[iclus];
        // float* _cluster_FEMC_Phi[iclus];
        // int* _cluster_FEMC_NTower[iclus];
        // float* _cluster_FEMC_trueID[iclus];
        for(Int_t iclus=0; iclus<_nclusters_FEMC; iclus++){
            if(verbosity>1) cout << "\tFEMC:  cluster " << iclus << "\twith E = " << _cluster_FEMC_E[iclus] << " GeV" << endl;
        }

        // ANCHOR FHCAL tower loop variables:
        // float* _tower_FHCAL_E[itwrH]
        // int* _tower_FHCAL_iEta[itwrH]
        // int* _tower_FHCAL_iPhi[itwrH]
        // float* _tower_FHCAL_trueID[itwrH]
        for(Int_t itwrH=0; itwrH<_nTowers_FHCAL; itwrH++){
            if(verbosity>1) cout << "\tFHCAL: tower " << itwrH << "\twith E = " << _tower_FHCAL_E[itwrH] << " GeV" << endl;
        }

        // ANCHOR FEMC cluster loop variables:
        // float* _tower_FEMC_E[itwrH]
        // int* _tower_FEMC_iEta[itwrH]
        // int* _tower_FEMC_iPhi[itwrH]
        // float* _tower_FEMC_trueID[itwrH]
        for(Int_t itwrE=0; itwrE<_nTowers_FEMC; itwrE++){
            if(verbosity>1) cout << "\tFEMC:  tower " << itwrE << "\twith E = " << _tower_FEMC_E[itwrE] << " GeV" << endl;
        }

        // ANCHOR MC particle loop variables:
        // float * _mcpart_ID[imc]
        // float * _mcpart_ID_parent[imc]
        // float * _mcpart_PDG[imc]
        // float * _mcpart_E [imc]
        // float * _mcpart_px[imc]
        // float * _mcpart_py[imc]
        // float * _mcpart_pz[imc]
        std::vector<float> jetf_truth_px;
        std::vector<float> jetf_truth_py;
        std::vector<float> jetf_truth_pz;
        std::vector<float> jetf_truth_E;
        for(Int_t imc=0; imc<_nMCPart; imc++){
            if(verbosity>1) cout << "\tMC:  particle " << imc << "\twith E = " << _mcpart_E[imc] << " GeV" << endl;
            //  create truth vector for jet finder
            jetf_truth_px.push_back(_mcpart_px[imc]);
            jetf_truth_py.push_back(_mcpart_py[imc]);
            jetf_truth_pz.push_back(_mcpart_pz[imc]);
            jetf_truth_E.push_back(_mcpart_E[imc]);
        }

        // ANCHOR JET FINDING
        // truth jets
        // auto jetsTrue = findJets(1.0, "anti-kt", jetf_truth_px, jetf_truth_py, jetf_truth_pz, jetf_truth_E);
        // if(verbosity>1)cout << "found " << std::get<1>(jetsTrue).size() << " true jets" << endl;        // printJets(std::get<1>(jetsTrue));

        // // track-based jets (rec)
        // auto jetsTrackRec = findJets(1.0, "anti-kt", jetf_track_px, jetf_track_py, jetf_track_pz, jetf_track_E);
        // if(verbosity>1)cout << "found " << std::get<1>(jetsTrackRec).size() << " rec track jets" << endl;        // printJets(std::get<1>(jetsTrackRec));

        // for (std::size_t i = 0; i < std::get<1>(jetsTrue).size(); i++) {
        //     // std::cout << "jet " << i << ": "<< std::get<1>(jetsTrue)[i].pt() << " " << std::get<1>(jetsTrue)[i].rap() << " " << std::get<1>(jetsTrue)[i].phi() << endl;
        //     for (std::size_t j = 0; j < std::get<1>(jetsTrackRec).size(); j++) {
        //         Double_t deltaRTrueRec = std::get<1>(jetsTrue)[i].delta_R(std::get<1>(jetsTrackRec)[j]);
        //         // cout << deltaRTrueRec << endl;
        //         if(deltaRTrueRec<0.5){
        //             if(verbosity>1) cout << "true jet " << i << " matched with rec jet " << j << " with dR = " << deltaRTrueRec << endl;
        //         }
        //     }
        // } // jet loop end
        resolutionhistos();

    } // event loop end
    resolutionhistosSave();
}
