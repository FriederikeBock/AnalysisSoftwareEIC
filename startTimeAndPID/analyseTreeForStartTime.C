//#include "/Users/syang/Tools/Macro/headers.h"
//#include "../common/plottingheader.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iomanip>

#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TVector3.h"

const Double_t mTiny = 1.e-6;
const Double_t mPtCutForT0 = 4;
//const Double_t mPtCutForT0 = 1000;

void analyseTreeForStartTime(
        TString studyname           = "ALLSILICON-FTTLSE2LC-ETTL-CTTLSE1-ACLGAD-TREXTOUT_epMB",
        TString singleFileName      = "",
        TString inputFileNameList   = "dataList/test.list",
        Bool_t isLANL               = kFALSE,
        Bool_t hasTiming            = kTRUE,
        Bool_t hasTrkSel            = kTRUE
        )
{
    TH1::SetDefaultSumw2(kTRUE);

    //gROOT->Reset();
    //gROOT->SetStyle("Plain");
    //StyleSettingsThesis();
    //TString dateForOutput           = ReturnDateStringForOutput();

    TString outputDir               = "rootfiles/";
    gSystem->Exec("mkdir -p "+outputDir);

    // **********************************************************************************************
    // ***************************** general variable setups ****************************************
    // **********************************************************************************************
    // random number generator
    TRandom3 r3;
    const double c                  = 29.9792458; // cm/ns
    const double sigmat             = 20e-3; // ns

    // regions: 0 = Backwards (electron going), 1 = central barrel, 2 = forward (hadron going)
    const Int_t nEtaRegs = 3;
    TString nameLayerLBL[nEtaRegs] = {"BACKWARD", "CENTRAL", "FORWARD"};
    TString nameLayerTTL[nEtaRegs] = {"ETTL", "CTTL", "FTTL"};

    Int_t maxLayer[nEtaRegs]       = {5, 6, 5};
    Int_t maxLayerTTL[nEtaRegs]    = {2, 2, 3};
    Int_t layerOffsetLBL[nEtaRegs] = {30, 10, 20};
    Int_t layerTTLBEC[nEtaRegs]    = {2, 1, 2};

    if (studyname.Contains("FTTLS2LC")){
        maxLayerTTL[2] = 2;
        layerTTLBEC[2] = 1;
    } else if( studyname.Contains("FTTLSE2LC")){
        maxLayerTTL[2] = 2;
        layerTTLBEC[2] = 2;
    } else if (studyname.Contains("FTTLSE1") ){
        maxLayerTTL[2] = 1;
        layerTTLBEC[2] = 1;
    }

    if (studyname.Contains("ETTLSE1")){
        maxLayerTTL[0] = 1;
        layerTTLBEC[0] = 1;
    }

    if (studyname.Contains("CTTLSE1")){
        maxLayerTTL[1] = 1;
        layerTTLBEC[1] = 1;
    }

    // **********************************************************************************************
    // ********************************** Set up histos *********************************************
    // **********************************************************************************************
    TH1D* hNEvents  = new TH1D("hNEvents", "nEvents", 1, 0, 1);
    TH1D* hNtrk    = new TH1D("hNtrk", "nTracks", 100, 0, 100);

    const   Int_t nSces = 2;
    TString sceName[nSces] = {"wScatE", "woScatE"};

    TH2D *hInitT0VsNtrk[nSces];
    TH2D *hInitT0VsNttlhit[nSces];
    TH2D *hNttlhitVsNtrk_wInitT0[nSces];
    TH2D *hNttlhitVsNtrk_wFinalT0[nSces];
    TH2D *hT0VsNtrk[nSces];
    TH2D *hT0VsNttlhit[nSces];
    TH2D *hInitDbetaVsP_pi[nSces];
    TH2D *hInitDbetaVsP_k[nSces];
    TH2D *hInitDbetaVsP_p[nSces];
    TH2D *hDbetaVsP_pi[nSces];
    TH2D *hDbetaVsP_k[nSces];
    TH2D *hDbetaVsP_p[nSces];
    TH2D *hBetaVsP_T0[nEtaRegs];
    TH2D *hDbetaVsP_wT0_piID[nEtaRegs];
    TH2D *hDbetaVsP_wT0_kID[nEtaRegs];
    TH2D *hDbetaVsP_wT0_pID[nEtaRegs];
    TH2D *hDbetaVsP_woT0_piID[nEtaRegs];
    TH2D *hDbetaVsP_woT0_kID[nEtaRegs];
    TH2D *hDbetaVsP_woT0_pID[nEtaRegs];
    TH2D *hDgenbetaVsGenP_piID[nEtaRegs];
    TH2D *hDgenbetaVsGenP_kID[nEtaRegs];
    TH2D *hDgenbetaVsGenP_pID[nEtaRegs];
    TH2D *hDbetaVsP_woT0_test_piID[nEtaRegs];
    TH2D *hDbetaVsP_woT0_test_kID[nEtaRegs];
    TH2D *hDbetaVsP_woT0_test_pID[nEtaRegs];
    for(int i=0; i<nSces; i++){
        hInitT0VsNtrk[i] = new TH2D(Form("hInitT0VsNtrk_%s", sceName[i].Data()), "hInitT0VsNtrk; N_{trk}; Initial T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);
        hInitT0VsNttlhit[i] = new TH2D(Form("hInitT0VsNttlhit_%s", sceName[i].Data()), "hInitT0VsNttlhit; N_{TTLHit}; Initial T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);

        hNttlhitVsNtrk_wInitT0[i]  = new TH2D(Form("hNttlhitVsNtrk_wInitT0_%s", sceName[i].Data()), "hNttlhitVsNtrk_wInitT0; N_{trk}; N_{TTLHit}", 100, 0, 100, 100, 0, 100);
        hNttlhitVsNtrk_wFinalT0[i] = new TH2D(Form("hNttlhitVsNtrk_wFinalT0_%s", sceName[i].Data()), "hNttlhitVsNtrk_wFinalT0; N_{trk}; N_{TTLHit}", 100, 0, 100, 100, 0, 100);

        hT0VsNtrk[i] = new TH2D(Form("hT0VsNtrk_%s", sceName[i].Data()), "hT0VsNtrk; N_{trk}; Final T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);
        hT0VsNttlhit[i] = new TH2D(Form("hT0VsNttlhit_%s", sceName[i].Data()), "hT0VsNttlhit; N_{TTLHit}; Final T_{0} (ns)", 100, 0, 100, 500, -0.5, 0.5);

        hInitDbetaVsP_pi[i] = new TH2D(Form("hInitDbetaVsP_pi_%s", sceName[i].Data()), "hInitDbetaVsP_pi; p (GeV); 1/#beta - 1/#beta_{#pi}", 300, 0, 15, 500, -0.5, 0.5);
        hInitDbetaVsP_k[i] = new TH2D(Form("hInitDbetaVsP_k_%s", sceName[i].Data()), "hInitDbetaVsP_k; p (GeV); 1/#beta - 1/#beta_{k}", 300, 0, 15, 500, -0.5, 0.5);
        hInitDbetaVsP_p[i] = new TH2D(Form("hInitDbetaVsP_p_%s", sceName[i].Data()), "hInitDbetaVsP_p; p (GeV); 1/#beta - 1/#beta_{p}", 300, 0, 15, 500, -0.5, 0.5);

        hDbetaVsP_pi[i] = new TH2D(Form("hDbetaVsP_pi_%s", sceName[i].Data()), "hDbetaVsP_pi; p (GeV); 1/#beta - 1/#beta_{#pi}", 300, 0, 15, 500, -0.5, 0.5);
        hDbetaVsP_k[i] = new TH2D(Form("hDbetaVsP_k_%s", sceName[i].Data()), "hDbetaVsP_k; p (GeV); 1/#beta - 1/#beta_{k}", 300, 0, 15, 500, -0.5, 0.5);
        hDbetaVsP_p[i] = new TH2D(Form("hDbetaVsP_p_%s", sceName[i].Data()), "hDbetaVsP_p; p (GeV); 1/#beta - 1/#beta_{p}", 300, 0, 15, 500, -0.5, 0.5);
    }

    for(Int_t eR=0; eR<nEtaRegs; eR++){
        hBetaVsP_T0[eR] = new TH2D(Form("hBetaVsP_T0_%s", nameLayerTTL[eR].Data()), "hBetaVsP_T0; p (GeV); 1/#beta", 300, 0, 15, 600, 0.9, 1.5);

        hDbetaVsP_wT0_piID[eR] = new TH2D(Form("hDbetaVsP_wT0_piID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_wT0_piID; p (GeV); 1/#beta - 1/#beta_{#pi}", 300, 0, 15, 500, -0.1, 0.1);
        hDbetaVsP_wT0_kID[eR] = new TH2D(Form("hDbetaVsP_wT0_kID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_wT0_kID; p (GeV); 1/#beta - 1/#beta_{k}", 300, 0, 15, 500, -0.1, 0.1);
        hDbetaVsP_wT0_pID[eR] = new TH2D(Form("hDbetaVsP_wT0_pID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_wT0_pID; p (GeV); 1/#beta - 1/#beta_{p}", 300, 0, 15, 500, -0.1, 0.1);

        hDbetaVsP_woT0_piID[eR] = new TH2D(Form("hDbetaVsP_woT0_piID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_woT0_piID; p (GeV); 1/#beta - 1/#beta_{#pi}", 300, 0, 15, 500, -0.1, 0.1);
        hDbetaVsP_woT0_kID[eR] = new TH2D(Form("hDbetaVsP_woT0_kID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_woT0_kID; p (GeV); 1/#beta - 1/#beta_{k}", 300, 0, 15, 500, -0.1, 0.1);
        hDbetaVsP_woT0_pID[eR] = new TH2D(Form("hDbetaVsP_woT0_pID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_woT0_pID; p (GeV); 1/#beta - 1/#beta_{p}", 300, 0, 15, 500, -0.1, 0.1);

        hDgenbetaVsGenP_piID[eR] = new TH2D(Form("hDgenbetaVsGenP_piID_%s", nameLayerTTL[eR].Data()), "hDgenbetaVsGenP_piID; p_{gen} (GeV); 1/#beta_{#pi}^{gen} - 1/#beta_{#pi}^{pathlength}", 300, 0, 15, 600, -0.3, 0.3);
        hDgenbetaVsGenP_kID[eR] = new TH2D(Form("hDgenbetaVsGenP_kID_%s", nameLayerTTL[eR].Data()), "hDgenbetaVsGenP_kID; p_{gen} (GeV); 1/#beta_{k}^{gen} - 1/#beta_{k}^{pathlength}", 300, 0, 15, 600, -0.3, 0.3);
        hDgenbetaVsGenP_pID[eR] = new TH2D(Form("hDgenbetaVsGenP_pID_%s", nameLayerTTL[eR].Data()), "hDgenbetaVsGenP_pID; p (GeV); p_{gen} (GeV); 1/#beta_{p}^{gen} - 1/#beta_{p}^{pathlength}", 300, 0, 15, 600, -0.3, 0.3);

        hDbetaVsP_woT0_test_piID[eR] = new TH2D(Form("hDbetaVsP_woT0_test_piID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_woT0_test_piID; p (GeV); 1/#beta - 1/#beta_{#pi}", 300, 0, 15, 500, -0.1, 0.1);
        hDbetaVsP_woT0_test_kID[eR] = new TH2D(Form("hDbetaVsP_woT0_test_kID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_woT0_test_kID; p (GeV); 1/#beta - 1/#beta_{k}", 300, 0, 15, 500, -0.1, 0.1);
        hDbetaVsP_woT0_test_pID[eR] = new TH2D(Form("hDbetaVsP_woT0_test_pID_%s", nameLayerTTL[eR].Data()), "hDbetaVsP_woT0_test_pID; p (GeV); 1/#beta - 1/#beta_{p}", 300, 0, 15, 500, -0.1, 0.1);
    }

    // **********************************************************************************************
    //************************** Read data **********************************************************
    // **********************************************************************************************
    TChain* t_tracks = new TChain("tracks", "tracks");
    if (singleFileName.CompareTo("") != 0){
        t_tracks->Add(singleFileName.Data());
    } else {
        ifstream in(inputFileNameList.Data());
        std::cout<<"File: "<<inputFileNameList.Data()<<std::endl;
        std::cout<<"reading input files "<<std::endl;
        // general number of triggers set
        while(!in.eof()){
            TString tempFilename = "";
            in >> tempFilename;
            if (tempFilename.CompareTo("") != 0 ){
                std::cout << "adding " << tempFilename.Data() << std::endl;
                t_tracks->Add(tempFilename.Data());
            }
        }      
    }

    if(t_tracks){
        Int_t eventID;
        Float_t trk_px,trk_py,trk_pz,trk_px_true,trk_py_true,trk_pz_true;
        Int_t pid;
        // recHit and trueHit arrays as follows: [etaRegion][layer][dimension] 
        // currently maximum number of layers = 6 (central barrel) organized inner to outer
        Float_t recHit[3][6][4];
        Float_t trueHit[3][6][4];
        // recHit and trueHit arrays as follows: [etaRegion][layer][dimension]
        // currently maximum number of layers = 3 (forward) organized inner to outer
        Float_t recHitTTL[3][3][4];
        Float_t trueHitTTL[3][3][4];

        t_tracks->SetBranchAddress("event", &eventID);
        t_tracks->SetBranchAddress("gflavor", &pid);
        t_tracks->SetBranchAddress("px", &trk_px);
        t_tracks->SetBranchAddress("py", &trk_py);
        t_tracks->SetBranchAddress("pz", &trk_pz);
        t_tracks->SetBranchAddress("gpx", &trk_px_true);
        t_tracks->SetBranchAddress("gpy", &trk_py_true);
        t_tracks->SetBranchAddress("gpz", &trk_pz_true);
        if (isLANL){
            for (Int_t l = 0; l < 5; l++){
                t_tracks->SetBranchAddress(Form("FST_%d_proj_x", l), &trueHit[2][l][0]);
                t_tracks->SetBranchAddress(Form("FST_%d_proj_y", l), &trueHit[2][l][1]);
                t_tracks->SetBranchAddress(Form("FST_%d_proj_z", l), &trueHit[2][l][2]);
                t_tracks->SetBranchAddress(Form("FST_%d_proj_t", l), &trueHit[2][l][3]); // pathlength
                t_tracks->SetBranchAddress(Form("FST_%d_x", l), &recHit[2][l][0]);
                t_tracks->SetBranchAddress(Form("FST_%d_y", l), &recHit[2][l][1]);
                t_tracks->SetBranchAddress(Form("FST_%d_z", l), &recHit[2][l][2]);
                t_tracks->SetBranchAddress(Form("FST_%d_t", l), &recHit[2][l][3]);
            }
        } 
        else {
            for (Int_t eR = 0; eR < nEtaRegs; eR++){
                for (Int_t l = 0; l < maxLayer[eR]; l++){
                    TString nameCurrBase = Form("LBLVTX_%s_%d", nameLayerLBL[eR].Data(), layerOffsetLBL[eR]+l);
                    std::cout << eR << "\t" << l << "\t" << nameCurrBase.Data() << std::endl;

                    t_tracks->SetBranchAddress(Form("%s_proj_x", nameCurrBase.Data()), &trueHit[eR][l][0]);
                    t_tracks->SetBranchAddress(Form("%s_proj_y", nameCurrBase.Data()), &trueHit[eR][l][1]);
                    t_tracks->SetBranchAddress(Form("%s_proj_z", nameCurrBase.Data()), &trueHit[eR][l][2]);
                    t_tracks->SetBranchAddress(Form("%s_proj_t", nameCurrBase.Data()), &trueHit[eR][l][3]); // pathlength
                    t_tracks->SetBranchAddress(Form("%s_x", nameCurrBase.Data()), &recHit[eR][l][0]);
                    t_tracks->SetBranchAddress(Form("%s_y", nameCurrBase.Data()), &recHit[eR][l][1]);
                    t_tracks->SetBranchAddress(Form("%s_z", nameCurrBase.Data()), &recHit[eR][l][2]);
                    t_tracks->SetBranchAddress(Form("%s_t", nameCurrBase.Data()), &recHit[eR][l][3]);      
                }
            }
        }

        if (hasTiming) {
            for (Int_t eR = 0; eR < nEtaRegs; eR++){
                for (Int_t l = 0;  l < maxLayerTTL[eR]; l++){
                    TString nameCurrBase = Form("%s_%d", nameLayerTTL[eR].Data(), l);
                    std::cout << eR << "\t" << l << "\t" << nameCurrBase.Data() << std::endl;

                    t_tracks->SetBranchAddress(Form("%s_proj_x", nameCurrBase.Data()), &trueHitTTL[eR][l][0]);
                    t_tracks->SetBranchAddress(Form("%s_proj_y", nameCurrBase.Data()), &trueHitTTL[eR][l][1]);
                    t_tracks->SetBranchAddress(Form("%s_proj_z", nameCurrBase.Data()), &trueHitTTL[eR][l][2]);
                    t_tracks->SetBranchAddress(Form("%s_proj_t", nameCurrBase.Data()), &trueHitTTL[eR][l][3]); // pathlength
                    t_tracks->SetBranchAddress(Form("%s_x", nameCurrBase.Data()), &recHitTTL[eR][l][0]);
                    t_tracks->SetBranchAddress(Form("%s_y", nameCurrBase.Data()), &recHitTTL[eR][l][1]);
                    t_tracks->SetBranchAddress(Form("%s_z", nameCurrBase.Data()), &recHitTTL[eR][l][2]);
                    t_tracks->SetBranchAddress(Form("%s_t", nameCurrBase.Data()), &recHitTTL[eR][l][3]);
                }
            }
        }

        std::vector<float> pathlength;
        std::vector<float> recMom, recTime, genMom, genTime;
        std::vector<int>   etaReg;
        std::vector<int>   trkID, uniqueTrkID;
        std::vector<bool>  isScatE; // scattering electron, assume they can be detected by e-going side EMC, which is directly idenfied by the geantID in this macro
        std::vector<int>   mcID;

        TParticlePDG *par = TDatabasePDG::Instance()->GetParticle(11);
        Double_t mElectron = par->Mass();

        par = TDatabasePDG::Instance()->GetParticle(211);
        Double_t mPion = par->Mass();

        par = TDatabasePDG::Instance()->GetParticle(321);
        Double_t mKaon = par->Mass();

        par = TDatabasePDG::Instance()->GetParticle(2212);
        Double_t mProton = par->Mass();

        const Int_t nIters = 4;
        float mT0[nIters];

        Int_t nentries = Int_t(t_tracks->GetEntries());
        cout<<"Total Entries: "<<nentries<<endl;

        Int_t currEventID = -1;
        Int_t nTrksCurr = 0;
        for (Long_t i = 0; i < nentries; ++i) {
            if (t_tracks->LoadTree(i) < 0) break;

            if(nentries<100){
                cout << "begin " << i << "th entry...." << endl;
            }
            else if(i % (nentries / 100) == 0){
                cout << "begin " << i << "th entry...." << endl;
            }

            t_tracks->GetEntry(i);

            if (eventID != currEventID) {
                if(i>0) { // deal with previous event
                    hNEvents->Fill(0.5);
                    hNtrk->Fill(nTrksCurr);

                    for(int it=0; it<nIters; it++) mT0[it] = -9999;

                    //cout<<"Before erase:"<<endl;
                    sort(uniqueTrkID.begin(), uniqueTrkID.end());
                    //for(size_t itrk=0; itrk<uniqueTrkID.size(); itrk++) cout<<uniqueTrkID[itrk]<<endl;
                    uniqueTrkID.erase(unique(uniqueTrkID.begin(), uniqueTrkID.end()), uniqueTrkID.end());
                    //cout<<"After erase:"<<endl;
                    //for(size_t itrk=0; itrk<uniqueTrkID.size(); itrk++) cout<<uniqueTrkID[itrk]<<endl;
                    //cout<<endl;

                    Int_t nScatEs = 0;
                    Int_t sceIdx = -1;
                    for(size_t ipar=0; ipar<recMom.size(); ipar++){
                        if(!isScatE[ipar]) continue;

                        float beta    = recMom[ipar] / sqrt(recMom[ipar]*recMom[ipar] + mElectron*mElectron);
                        float tflight = pathlength[ipar] / (beta * c);

                        if(nScatEs == 0) mT0[0]  = recTime[ipar] - tflight;
                        else             mT0[0] += recTime[ipar] - tflight;

                        nScatEs++;
                    }

                    if(nScatEs > 0){ // evaluate initial T0 with the scattering electron
                        sceIdx = 0;
                        mT0[0] = mT0[0] / nScatEs;
                    }
                    else{ // All particles are assumed to be pion for the events without scattering electron detection to evaluate initial T0
                        for(size_t ipar=0; ipar<recMom.size(); ipar++){
                            float beta    = recMom[ipar] / sqrt(recMom[ipar]*recMom[ipar] + mPion*mPion);
                            float tflight = pathlength[ipar] / (beta * c);

                            if(ipar == 0) mT0[0]  = recTime[ipar] - tflight;
                            else          mT0[0] += recTime[ipar] - tflight;
                        }

                        if(recTime.size() > 0){
                            sceIdx = 1;
                            mT0[0] = mT0[0] / (recTime.size());
                        }
                    }

                    if(sceIdx > -1){
                        hInitT0VsNtrk[sceIdx]->Fill(nTrksCurr, mT0[0]);
                        hInitT0VsNttlhit[sceIdx]->Fill(trkID.size(), mT0[0]);
                    }

                    // do iteration to identify particles to improve T0 determination
                    for(Int_t it=1; it<nIters; it++){
                        Int_t  nUsedTTLHits = 0;

                        if(fabs(mT0[it-1] + 9999) < mTiny) break;

                        for(size_t ipar=0; ipar<recMom.size(); ipar++){
                            Bool_t isPID = kFALSE;

                            float tflight   = recTime[ipar] - mT0[it-1];
                            float beta      = pathlength[ipar] / tflight / c;
                            float beta_woT0 = pathlength[ipar] / recTime[ipar] / c;

                            float beta_e  = recMom[ipar] / sqrt(recMom[ipar]*recMom[ipar] + mElectron*mElectron);
                            float beta_pi = recMom[ipar] / sqrt(recMom[ipar]*recMom[ipar] + mPion*mPion);
                            float beta_k  = recMom[ipar] / sqrt(recMom[ipar]*recMom[ipar] + mKaon*mKaon);
                            float beta_p  = recMom[ipar] / sqrt(recMom[ipar]*recMom[ipar] + mProton*mProton);

                            float dBeta_pi = 1/beta - 1/beta_pi;
                            float dBeta_k  = 1/beta - 1/beta_k;
                            float dBeta_p  = 1/beta - 1/beta_p;

                            float genBeta    = pathlength[ipar] / genTime[ipar] / c;
                            float genBeta_pi = genMom[ipar] / sqrt(genMom[ipar]*genMom[ipar] + mPion*mPion);
                            float genBeta_k  = genMom[ipar] / sqrt(genMom[ipar]*genMom[ipar] + mKaon*mKaon);
                            float genBeta_p  = genMom[ipar] / sqrt(genMom[ipar]*genMom[ipar] + mProton*mProton);

                            if(!isScatE[ipar] && it==1){
                                hInitDbetaVsP_pi[sceIdx]->Fill(recMom[ipar], dBeta_pi);
                                hInitDbetaVsP_k[sceIdx]->Fill(recMom[ipar], dBeta_k);
                                hInitDbetaVsP_p[sceIdx]->Fill(recMom[ipar], dBeta_p);

                                if(fabs(mcID[ipar]) == 211){ // pion
                                    hDgenbetaVsGenP_piID[etaReg[ipar]]->Fill(genMom[ipar], 1/genBeta - 1/genBeta_pi);
                                    hDbetaVsP_woT0_test_piID[etaReg[ipar]]->Fill(recMom[ipar], 1/beta_woT0 - 1/beta_pi);
                                }

                                if(fabs(mcID[ipar]) == 321){ // kaon
                                    hDgenbetaVsGenP_kID[etaReg[ipar]]->Fill(genMom[ipar], 1/genBeta - 1/genBeta_k);
                                    hDbetaVsP_woT0_test_kID[etaReg[ipar]]->Fill(recMom[ipar], 1/beta_woT0 - 1/beta_k);
                                }

                                if(fabs(mcID[ipar]) == 2212){ // prton
                                    hDgenbetaVsGenP_pID[etaReg[ipar]]->Fill(genMom[ipar], 1/genBeta - 1/genBeta_p);
                                    hDbetaVsP_woT0_test_pID[etaReg[ipar]]->Fill(recMom[ipar], 1/beta_woT0 - 1/beta_p);
                                }
                            }

                            if(isScatE[ipar]){
                                isPID   = kTRUE;
                                tflight = pathlength[ipar] / (beta_e * c);

                                nUsedTTLHits++;
                            }
                            else if(recMom[ipar] < mPtCutForT0 && dBeta_pi < 0.005 && dBeta_pi > -0.01){
                                isPID   = kTRUE;
                                tflight = pathlength[ipar] / (beta_pi * c);

                                nUsedTTLHits++;
                            }
                            else if(recMom[ipar] < mPtCutForT0 && fabs(dBeta_k) < 0.005){
                                isPID   = kTRUE;
                                tflight = pathlength[ipar] / (beta_k * c);

                                nUsedTTLHits++;
                            }
                            else if(recMom[ipar] < mPtCutForT0 && fabs(dBeta_p) < 0.005){
                                isPID   = kTRUE;
                                tflight = pathlength[ipar] / (beta_p * c);

                                nUsedTTLHits++;
                            }

                            if(isPID){
                                if(fabs(mT0[it] + 9999) < mTiny) mT0[it] = recTime[ipar] - tflight;
                                else mT0[it] += recTime[ipar] - tflight;
                            }
                        }
                        //cout<<nUsedTTLHits<<endl;

                        if(nUsedTTLHits > 0){
                            mT0[it] = mT0[it]/nUsedTTLHits;
                        }

                        //for(int it=0; it<nIters; it++){
                        //    cout<<"Iter = " << it << "   : T0 = "<< mT0[it]<<endl; 
                        //}
                        //cout<<endl;
                    }

                    float mFinalT0 = -9999;
                    if(nScatEs > 0){
                        for(Int_t it=0; it<nIters; it++){ // for the events with the scattering electron detected, the initial T0 could be used as the final T0
                            if(mT0[it] > -9999) mFinalT0 = mT0[it];
                        }
                    }
                    else{
                        for(Int_t it=1; it<nIters; it++){
                            if(mT0[it] > -9999) mFinalT0 = mT0[it];
                        }
                    }

                    //if(mT0[0] > -9999 && fabs(mFinalT0 + 9999)<mTiny) cout<<"T0 calibration failed"<<endl;
                    if(mT0[0] > -9999){
                        hNttlhitVsNtrk_wInitT0[sceIdx]->Fill(nTrksCurr, trkID.size());
                        if(mFinalT0 > -9999) hNttlhitVsNtrk_wFinalT0[sceIdx]->Fill(nTrksCurr, trkID.size());
                    }

                    if(mFinalT0 > -9999){
                        hT0VsNtrk[sceIdx]->Fill(nTrksCurr, mFinalT0);
                        hT0VsNttlhit[sceIdx]->Fill(trkID.size(), mFinalT0);

                        int   currTrkID = -1;
                        float totalPathLength = 0, totalTFlight = 0, totalTFlightWoT0 = 0;
                        for(size_t ipar=0; ipar<recMom.size(); ipar++){
                            if(trkID[ipar] != currTrkID){
                                if(ipar > 0){ // deal with previous track
                                    float beta = totalPathLength / totalTFlight / c;
                                    hBetaVsP_T0[etaReg[ipar-1]]->Fill(recMom[ipar-1], 1/beta);

                                    float beta_pi = recMom[ipar-1] / sqrt(recMom[ipar-1]*recMom[ipar-1] + mPion*mPion);
                                    float beta_k  = recMom[ipar-1] / sqrt(recMom[ipar-1]*recMom[ipar-1] + mKaon*mKaon);
                                    float beta_p  = recMom[ipar-1] / sqrt(recMom[ipar-1]*recMom[ipar-1] + mProton*mProton);

                                    float dBeta_pi = 1/beta - 1/beta_pi;
                                    float dBeta_k  = 1/beta - 1/beta_k;
                                    float dBeta_p  = 1/beta - 1/beta_p;

                                    hDbetaVsP_pi[sceIdx]->Fill(recMom[ipar-1], dBeta_pi);
                                    hDbetaVsP_k[sceIdx]->Fill(recMom[ipar-1], dBeta_k);
                                    hDbetaVsP_p[sceIdx]->Fill(recMom[ipar-1], dBeta_p);

                                    float beta_woT0     = totalPathLength / totalTFlightWoT0 / c;
                                    float dBeta_woT0_pi = 1/beta_woT0 - 1/beta_pi;
                                    float dBeta_woT0_k  = 1/beta_woT0 - 1/beta_k;
                                    float dBeta_woT0_p  = 1/beta_woT0 - 1/beta_p;

                                    if(fabs(mcID[ipar-1]) == 211){ // pion
                                        hDbetaVsP_wT0_piID[etaReg[ipar-1]]->Fill(recMom[ipar-1], dBeta_pi);
                                        hDbetaVsP_woT0_piID[etaReg[ipar-1]]->Fill(recMom[ipar-1], dBeta_woT0_pi);
                                    }

                                    if(fabs(mcID[ipar-1]) == 321){ // kaon
                                        hDbetaVsP_wT0_kID[etaReg[ipar-1]]->Fill(recMom[ipar-1], dBeta_k);
                                        hDbetaVsP_woT0_kID[etaReg[ipar-1]]->Fill(recMom[ipar-1], dBeta_woT0_k);
                                    }

                                    if(fabs(mcID[ipar-1]) == 2212){ // prton
                                        hDbetaVsP_wT0_pID[etaReg[ipar-1]]->Fill(recMom[ipar-1], dBeta_p);
                                        hDbetaVsP_woT0_pID[etaReg[ipar-1]]->Fill(recMom[ipar-1], dBeta_woT0_p);
                                    }
                                }

                                currTrkID = trkID[ipar];
                                totalPathLength = 0; totalTFlight = 0, totalTFlightWoT0 = 0;

                                totalPathLength  += pathlength[ipar];
                                totalTFlight     += recTime[ipar] - mFinalT0;
                                totalTFlightWoT0 += recTime[ipar];
                            }
                            else{
                                totalPathLength  += pathlength[ipar];
                                totalTFlight     += recTime[ipar] - mFinalT0;
                                totalTFlightWoT0 += recTime[ipar];
                            }
                        }

                        int   pIdx = recMom.size()-1;
                        float beta = totalPathLength / totalTFlight / c;
                        hBetaVsP_T0[etaReg[pIdx]]->Fill(recMom[pIdx], 1/beta);

                        float beta_pi = recMom[pIdx] / sqrt(recMom[pIdx]*recMom[pIdx] + mPion*mPion);
                        float beta_k  = recMom[pIdx] / sqrt(recMom[pIdx]*recMom[pIdx] + mKaon*mKaon);
                        float beta_p  = recMom[pIdx] / sqrt(recMom[pIdx]*recMom[pIdx] + mProton*mProton);

                        float dBeta_pi = 1/beta - 1/beta_pi;
                        float dBeta_k  = 1/beta - 1/beta_k;
                        float dBeta_p  = 1/beta - 1/beta_p;

                        hDbetaVsP_pi[sceIdx]->Fill(recMom[pIdx], dBeta_pi);
                        hDbetaVsP_k[sceIdx]->Fill(recMom[pIdx], dBeta_k);
                        hDbetaVsP_p[sceIdx]->Fill(recMom[pIdx], dBeta_p);

                        float beta_woT0     = totalPathLength / totalTFlightWoT0 / c;
                        float dBeta_woT0_pi = 1/beta_woT0 - 1/beta_pi;
                        float dBeta_woT0_k  = 1/beta_woT0 - 1/beta_k;
                        float dBeta_woT0_p  = 1/beta_woT0 - 1/beta_p;

                        if(fabs(mcID[pIdx]) == 211){ // pion
                            hDbetaVsP_wT0_piID[etaReg[pIdx]]->Fill(recMom[pIdx], dBeta_pi);
                            hDbetaVsP_woT0_piID[etaReg[pIdx]]->Fill(recMom[pIdx], dBeta_woT0_pi);
                        }

                        if(fabs(mcID[pIdx]) == 321){ // kaon
                            hDbetaVsP_wT0_kID[etaReg[pIdx]]->Fill(recMom[pIdx], dBeta_k);
                            hDbetaVsP_woT0_kID[etaReg[pIdx]]->Fill(recMom[pIdx], dBeta_woT0_k);
                        }

                        if(fabs(mcID[pIdx]) == 2212){ // prton
                            hDbetaVsP_wT0_pID[etaReg[pIdx]]->Fill(recMom[pIdx], dBeta_p);
                            hDbetaVsP_woT0_pID[etaReg[pIdx]]->Fill(recMom[pIdx], dBeta_woT0_p);
                        }
                    }
                }

                pathlength.clear();
                recMom.clear();
                recTime.clear();
                mcID.clear();
                genMom.clear();
                genTime.clear();
                etaReg.clear();
                trkID.clear(); uniqueTrkID.clear();
                isScatE.clear();

                currEventID = eventID;
                nTrksCurr = 1;
            }
            else {
                nTrksCurr++; 
            }

            if(!std::isnan(trk_px) && hasTiming){
                TVector3 rec_mom(trk_px, trk_py, trk_pz);
                TVector3 gen_mom(trk_px_true, trk_py_true, trk_pz_true);

                for(Int_t eR=0; eR<nEtaRegs; eR++){
                    for(Int_t l=0; l<maxLayerTTL[eR]; l++){
                        //for(Int_t l=0; l<layerTTLBEC[eR]; l++) // only use the TTL layers before EMC
                        if(recHitTTL[eR][l][3] > -9999 && trueHitTTL[eR][l][3]>0){ // has time information and pathlength>0
                            if(hasTrkSel){
                                TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(pid);
                                Double_t mass = part->Mass();

                                float expBeta = trueHitTTL[eR][l][3] / recHitTTL[eR][l][3] / c; // pathlength / trueFlightTime
                                float genBeta = gen_mom.Mag() / sqrt(gen_mom.Mag()*gen_mom.Mag() + mass*mass);
                                if(fabs(1/genBeta - 1/expBeta) >= 1.5e-2) continue; // reject wrong trajectory, note this is just a temporary selection, potentioally bias the results !!! 
                            }

                            pathlength.push_back(trueHitTTL[eR][l][3]);
                            recMom.push_back(rec_mom.Mag());
                            recTime.push_back(r3.Gaus(recHitTTL[eR][l][3], sigmat));
                            mcID.push_back(pid);
                            genMom.push_back(gen_mom.Mag());
                            genTime.push_back(recHitTTL[eR][l][3]);

                            etaReg.push_back(eR);
                            trkID.push_back(i); uniqueTrkID.push_back(i);


                            if(eR == 0 && fabs(pid) == 11) isScatE.push_back(true);
                            else                           isScatE.push_back(false);
                        }
                    }
                }
            }
        }

        // deal with the last event
        hNEvents->Fill(0.5);
        hNtrk->Fill(nTrksCurr);
        // suppose to evaluate the T0 for the last event here, but ....
    }

    TString nameOutput;
    if(inputFileNameList.Contains("test.list")){
        if(hasTrkSel) nameOutput = Form("%s/startTimePerformance_test.root", outputDir.Data());
        else          nameOutput = Form("%s/startTimePerformance_woTrkSel_test.root", outputDir.Data());
    }
    else if(mPtCutForT0 > 4){
        if(hasTrkSel) nameOutput = Form("%s/startTimePerformance_wHiPtParticles_%s.root",outputDir.Data(), studyname.Data());
        else          nameOutput = Form("%s/startTimePerformance_wHiPtParticles_woTrkSel_%s.root",outputDir.Data(), studyname.Data());
    }
    else{
        if(hasTrkSel) nameOutput = Form("%s/startTimePerformance_%s.root",outputDir.Data(), studyname.Data());
        else          nameOutput = Form("%s/startTimePerformance_woTrkSel_%s.root",outputDir.Data(), studyname.Data());
    }
    TFile* fOutput          = new TFile(nameOutput,"RECREATE");
    hNEvents->Write();
    hNtrk->Write();
    for(Int_t isce=0; isce<nSces; isce++){
        hInitT0VsNtrk[isce]->Write();
        hInitT0VsNttlhit[isce]->Write();
        hNttlhitVsNtrk_wInitT0[isce]->Write();
        hNttlhitVsNtrk_wFinalT0[isce]->Write();
        hT0VsNtrk[isce]->Write();
        hT0VsNttlhit[isce]->Write();
        hInitDbetaVsP_pi[isce]->Write();
        hInitDbetaVsP_k[isce]->Write();
        hInitDbetaVsP_p[isce]->Write();
        hDbetaVsP_pi[isce]->Write();
        hDbetaVsP_k[isce]->Write();
        hDbetaVsP_p[isce]->Write();
    }
    for(Int_t eR=0; eR<nEtaRegs; eR++){
        hBetaVsP_T0[eR]->Write();

        hDbetaVsP_wT0_piID[eR]->Write();
        hDbetaVsP_wT0_kID[eR]->Write();
        hDbetaVsP_wT0_pID[eR]->Write();

        hDbetaVsP_woT0_piID[eR]->Write();
        hDbetaVsP_woT0_kID[eR]->Write();
        hDbetaVsP_woT0_pID[eR]->Write();

        hDgenbetaVsGenP_piID[eR]->Write();
        hDgenbetaVsGenP_kID[eR]->Write();
        hDgenbetaVsGenP_pID[eR]->Write();

        hDbetaVsP_woT0_test_piID[eR]->Write();
        hDbetaVsP_woT0_test_kID[eR]->Write();
        hDbetaVsP_woT0_test_pID[eR]->Write();
    }

    fOutput->Close();
    delete fOutput;
}
