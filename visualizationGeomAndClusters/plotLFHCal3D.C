#include <algorithm>
#include "../common/plottingheader.h"
#include "../treeAnalysis/treeProcessing.h"
#include "../treeAnalysis/clusterizer.cxx"
// #include "../treeAnalysis/caloheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))
#include "TColor.h"
#include <TChain.h>
#include <TVector3.h>

typedef struct {
  float tower_Eta;
  float tower_Phi;
  float tower_x;
  float tower_y;
  float tower_z;
  float tower_R;
  float tower_theta;
} geomStrct;

void plotLFHCal3D(
    TString nameFile            = "/home/ewa/alice/EIC/DST_HFandJets_pythia8_ep-18x275-q2-100_000_0124000_02000_g4event_eval.root",
    TString inFileGeometry      = "/home/ewa/alice/EIC/geometry_ideal_BECALproj_EEMC_LFHCAL_corrected.root",
    TString suffix              = "pdf",
    TString collisionsSys       = "PythiaMB",
    double energy               = -1.0,
    int plotEvent               = 7,
    int clusterizerEnum         = 0,
    TString addName             = ""           
){
    // dimensions of dummy histo
    int towersx         = 3*106;          // 105 towers for LFHCAl, 3*105 for FEMC
    int towersy         = 3*106;          // 105 towers for LFHCAl, 3*105 for FEMC
    int towersz         = 8*7 + 15;       // 7 towers for LFHCAL, 1 for FEMC 
    // dimensions of LFHCAL
    int towersxyLF      = 105;
    int towerszLF       = 7;

    // histogram spacing in cm
    double dX           = 1.61;
    double dY           = 1.61;
    double dZ           = 2.5;


    TString labelDet  = Form("LFHCal%s", addName.Data());

    float minECutOff_LFHCAL = aggE[8]/2.;
    float seedEC_LFHCAL     = seedE[8];
    float aggEC_LFHCAL      = aggE[8];

    float minECutOff_FEMC   = aggE[1]/2.;
    float seedEC_FEMC       = seedE[1];
    float aggEC_FEMC        = aggE[1];

    int verbosity     = 1;

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    StyleSettingsThesis();

    Color_t colorCluster[50] = {kRed, kGreen+2, kBlue, kOrange, kCyan+1, kMagenta+1, kViolet-4, kTeal, kSpring-1, kYellow,
                            kRed+3,  kGreen+3, kBlue+3, kOrange+2, kCyan+3, kMagenta+3, kViolet+2, kTeal-1, kSpring+2, kYellow+2,
                            kRed-4, kGreen-4, kBlue-4, kOrange-5, kCyan-5, kMagenta-3, kViolet-3, kTeal-6, kSpring-6, kYellow-6,
                            kRed-8, kGreen-6, kBlue-7, kOrange+7, kCyan-6, kMagenta-6, kViolet-7, kTeal+3, kSpring+9, kYellow-5,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
    TString labelEnergy = "";
    if (collisionsSys.CompareTo("PythiaMB") == 0){
        labelEnergy   = "e-p: 18#times 275 GeV";
    } else if (collisionsSys.CompareTo("PythiaQ2_100") == 0){
        labelEnergy   = "e-p: 18#times 275 GeV, Q^{2} > 100 GeV";
    } else if (collisionsSys.CompareTo("Satre18x108") == 0){
        labelEnergy   = "e-p: 18#times 108 GeV, Satre";
    } else if (collisionsSys.CompareTo("SingleElectron") == 0){
        labelEnergy   = Form("e^{-}, %0.1f GeV", energy);
    } else if (collisionsSys.CompareTo("SingleE") == 0){
        labelEnergy   = "e^{-} simulation";
    } else if (collisionsSys.CompareTo("SingleKaon") == 0){
        labelEnergy   = Form("K^{-}, %0.1f GeV", energy);
    } else if (collisionsSys.CompareTo("SinglePion") == 0){
        labelEnergy   = Form("#pi^{-}, %0.1f GeV", energy);
    } else if (collisionsSys.CompareTo("SingleProton") == 0){
        labelEnergy   = Form("p, %0.1f GeV", energy);
    } else if (collisionsSys.CompareTo("SingleNeutron") == 0){
        labelEnergy   = Form("n, %0.1f GeV", energy);
    }
  
    TString dateForOutput                         = ReturnDateStringForOutput();
    TString outputDir                             = Form("plots/%s_%s_%s", dateForOutput.Data(), collisionsSys.Data(), labelDet.Data());
    if (collisionsSys.Contains("Single"))
        outputDir = Form("%s_%.0f", outputDir.Data(), energy);
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir+"/3D");
    gSystem->Exec("mkdir -p "+outputDir+"/3DwClusters");

    double  textSizeSinglePad                    = 0.05;
    double  textSizeLabelsPixel                  = 35;
    double  textSizeLabelsRel                    = 58./1300;

    int granularity = 1;
    TH3F*   h_XYZMapEvt                 = NULL;
    TH3F*   h_XYZMapEvt_Hole            = NULL;
    TH3F*   h_XYZMapEvt_Ecut            = NULL;
    TH3F*   h_XYZMapEvt_Cl[50]          = {NULL};
    TH3F*   h_XYZMapEvt_MCPart[50]      = {NULL};

    // reading file
    // load event tree
    TChain *tt_event  = new TChain("event_tree");
    tt_event->Add(nameFile.Data());
    if( !tt_event ) { std::cout << "tree not found... returning!" << std::endl; return; }
  
    // load geometry 
    tt_geometry   = (TTree*) (new TFile(inFileGeometry.Data(),"READ"))->Get("geometry_tree");
    if( !tt_geometry ) { std:: cout << "geometry tree not found... returning! " << std:: endl; return; }
  
    // load all branches (see treeProcessing header)
    SetBranchAddressesTree(tt_event);
    SetBranchAddressesGeometryTree(tt_geometry);
    SetGeometryIndices();

    // define the histos
    h_XYZMapEvt              = new TH3F("h_XYZMapEvt","",towersz, 312.5, 490, towersx, -262.1, 262.1, towersy, -262.1, 262.1);
    h_XYZMapEvt->Sumw2();
    h_XYZMapEvt_Hole         = new TH3F("h_XYZMapEvt_Hole","",towersz, 312.5, 490, towersx, -262.1, 262.1, towersy, -262.1, 262.1);
    h_XYZMapEvt_Hole->Sumw2();
    h_XYZMapEvt_Ecut         = new TH3F("h_XYZMapEvt_Ecut","",towersz, 312.5, 490, towersx, -262.1, 262.1, towersy, -262.1, 262.1);
    h_XYZMapEvt_Ecut->Sumw2();
    for(int i=0; i<50; i++){
        h_XYZMapEvt_Cl[i]       = new TH3F(Form("h_XYZMapEvt_Cl%i",i),"",towersz, 312.5, 490, towersx, -262.1, 262.1, towersy, -262.1, 262.1);
        h_XYZMapEvt_Cl[i]->Sumw2();
        h_XYZMapEvt_MCPart[i]   = new TH3F(Form("h_XYZMapEvt_MCPart%i",i),"",towersz, 312.5, 490, towersx, -262.1, 262.1, towersy, -262.1, 262.1);
        h_XYZMapEvt_MCPart[i]->Sumw2();
    }

    // load geometry of LFHCAL - caloID 8
    tt_geometry->GetEvent( calogeomindex[8] );
    int maxiPhi = 0;
    int maxiEta = 0;
    // vector geometry structure
    geomStrct   arr_LFHCALGeoStrt[towersxyLF][towersxyLF][towerszLF];  // array[iEta-2][iPhi-2]
    for(int i=0; i<(int)_calogeom_towers_N; i++){
        int tmpEta  = _calogeom_towers_iEta[i]-2;
        int tmpPhi  = _calogeom_towers_iPhi[i]-2;
        int tmpIL   = _calogeom_towers_iL[i];
        if (_calogeom_towers_iEta[i] > maxiEta) maxiEta = _calogeom_towers_iEta[i];
        if (_calogeom_towers_iPhi[i] > maxiPhi) maxiPhi = _calogeom_towers_iPhi[i];
        arr_LFHCALGeoStrt[tmpEta][tmpPhi][tmpIL].tower_Eta   = _calogeom_towers_Eta[i];
        arr_LFHCALGeoStrt[tmpEta][tmpPhi][tmpIL].tower_Phi   = _calogeom_towers_Phi[i];
        arr_LFHCALGeoStrt[tmpEta][tmpPhi][tmpIL].tower_x     = _calogeom_towers_x[i];
        arr_LFHCALGeoStrt[tmpEta][tmpPhi][tmpIL].tower_y     = _calogeom_towers_y[i];
        arr_LFHCALGeoStrt[tmpEta][tmpPhi][tmpIL].tower_z     = _calogeom_towers_z[i];
        for(int iX=0; iX<3; iX++){
            for(int iY=0; iY<3; iY++){
                for(int iZ=0; iZ<8; iZ++){
                    h_XYZMapEvt_Hole->SetBinContent( 8*(tmpIL)+iZ+16, 3*tmpEta+iX+3, 3*tmpPhi+iY+3, 1);
                }
            }
        }
    }
    std::cout << maxiEta << "\t" << maxiPhi << std::endl;

    // load geometry of LFHCAL - caloID 8
    tt_geometry->GetEvent( calogeomindex[1] );

    geomStrct   arr_FEMCGeoStrt[towersx][towersy];  
    for(int i=0; i<(int)_calogeom_towers_N; i++){
        int tmpEta  = _calogeom_towers_iEta[i]-2;
        int tmpPhi  = _calogeom_towers_iPhi[i]-2;
        arr_FEMCGeoStrt[tmpEta][tmpPhi].tower_Eta   = _calogeom_towers_Eta[i];
        arr_FEMCGeoStrt[tmpEta][tmpPhi].tower_Phi   = _calogeom_towers_Phi[i];
        arr_FEMCGeoStrt[tmpEta][tmpPhi].tower_x     = _calogeom_towers_x[i];
        arr_FEMCGeoStrt[tmpEta][tmpPhi].tower_y     = _calogeom_towers_y[i];
        arr_FEMCGeoStrt[tmpEta][tmpPhi].tower_z     = _calogeom_towers_z[i];
        for(int iZ=0; iZ<15; iZ++){
            h_XYZMapEvt_Hole->SetBinContent(iZ+1, tmpEta+62, tmpPhi+62, 1);
        }
    }

    for(int z = 1; z < h_XYZMapEvt_Hole->GetNbinsX()+1; z++){
        for(int x = 1; x < h_XYZMapEvt_Hole->GetNbinsY()+1; x++){
            for(int y = 1; y < h_XYZMapEvt_Hole->GetNbinsZ()+1; y++){
                if( z < 16 && ((x < 130 || x > 190) || (y < 130 || y > 190))){
                    h_XYZMapEvt_Hole->SetBinContent(z,x,y, 0); 
                } else {
                    if( h_XYZMapEvt_Hole->GetBinContent(z,x,y) > 0)
                        h_XYZMapEvt_Hole->SetBinContent(z,x,y, 0);
                    else 
                        h_XYZMapEvt_Hole->SetBinContent(z,x,y, 0.1);
                }
            }
        }
    }

    float   sumedE_LFHCAL       = 0;
    float   sumedE_FEMC         = 0;

    Int_t nentries  = Int_t(tt_event->GetEntries());
    for(Int_t i=0; i<nentries; i++){

        if( tt_event->LoadTree(i) < 0 ) break;
        if( i!=plotEvent ) continue;
        tt_event->GetEntry(i);

        for(Int_t imc=0; (Int_t)imc<_nMCPart; imc++){
          TVector3 truevec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
          _mcpart_Eta[imc]=truevec.Eta();
          _mcpart_Phi[imc]=truevec.Phi();
        }

        if( _nTowers_LFHCAL == 0 && _nTowers_FEMC == 0 ) return;

        // LFHCAL hits
        for(int itwr=0;(int)itwr<_nTowers_LFHCAL;itwr++){
            sumedE_LFHCAL += _tower_LFHCAL_E[itwr];
            if (_tower_LFHCAL_E[itwr] > 1e-3 && verbosity > 1) cout << _tower_LFHCAL_iEta[itwr] << "\t" << _tower_LFHCAL_iPhi[itwr] << "\t" << _tower_LFHCAL_iL[itwr] << "\t"<< _tower_LFHCAL_E[itwr] << endl;
            int tmpX = arr_LFHCALGeoStrt[_tower_LFHCAL_iEta[itwr]][_tower_LFHCAL_iPhi[itwr]][_tower_LFHCAL_iL[itwr]].tower_x;
            int tmpY = arr_LFHCALGeoStrt[_tower_LFHCAL_iEta[itwr]][_tower_LFHCAL_iPhi[itwr]][_tower_LFHCAL_iL[itwr]].tower_y;
            int tmpZ = arr_LFHCALGeoStrt[_tower_LFHCAL_iEta[itwr]][_tower_LFHCAL_iPhi[itwr]][_tower_LFHCAL_iL[itwr]].tower_z;
            h_XYZMapEvt->Fill(tmpZ + 3*dZ, tmpX + dX, tmpY + dY, _tower_LFHCAL_E[itwr]);
            if(_tower_LFHCAL_E[itwr]>minECutOff_LFHCAL){
                h_XYZMapEvt_Ecut->Fill(tmpZ + 3*dZ, tmpX + dX, tmpY + dY, _tower_LFHCAL_E[itwr]);
            }
        }

        // FEMC hits
        for(int itwr=0;(int)itwr<_nTowers_FEMC;itwr++){
            sumedE_FEMC += _tower_FEMC_E[itwr];
            if (_tower_FEMC_E[itwr] > 1e-3 && verbosity > 1) cout << _tower_FEMC_iEta[itwr] << "\t" << _tower_FEMC_iPhi[itwr] << "\t"<< _tower_FEMC_E[itwr] << endl;
            int tmpX = arr_FEMCGeoStrt[_tower_FEMC_iEta[itwr]][_tower_FEMC_iPhi[itwr]].tower_x;
            int tmpY = arr_FEMCGeoStrt[_tower_FEMC_iEta[itwr]][_tower_FEMC_iPhi[itwr]].tower_y;
            int tmpZ = arr_FEMCGeoStrt[_tower_FEMC_iEta[itwr]][_tower_FEMC_iPhi[itwr]].tower_z;
            h_XYZMapEvt->Fill(tmpZ + 7*dZ, tmpX, tmpY, _tower_FEMC_E[itwr]);
            if(_tower_FEMC_E[itwr]>minECutOff_FEMC){
                h_XYZMapEvt_Ecut->Fill(tmpZ + 7*dZ, tmpX, tmpY, _tower_FEMC_E[itwr]);
            }
        }
    }

    int nclusters = 0;

    // ************************************************************************************************************************************************************
    // LFHCAL clusterizer
    
    // vector of towers used in the clusterizer
    std::vector<towersStrct> input_towers_LFHCAL;
    std::vector<towersStrct> input_towersForMC_LFHCAL;
    // vector of towers within the currently found cluster
    std::vector<towersStrct> cluster_towers_LFHCAL;
    std::vector<clustersStrct> rec_clusters_LFHCAL;
    std::vector<float> recclusterEnergies_LFHCAL;
    std::vector<TString> mcPartEnergies;

    rec_clusters_LFHCAL.clear();
    recclusterEnergies_LFHCAL.clear();
    mcPartEnergies.clear();

    int caloID = 8;     // clusterization for LFHCAL

    // fill vector with towers for clusterization above aggregation threshold
    for(int itow=0; itow<(int)_nTowers_LFHCAL; itow++){
        if(_tower_LFHCAL_E[itow]>aggEC_LFHCAL){
            towersStrct tempstructT;
            tempstructT.tower_E         = _tower_LFHCAL_E[itow];
            tempstructT.tower_iEta      = _tower_LFHCAL_iEta[itow];
            tempstructT.tower_iPhi      = _tower_LFHCAL_iPhi[itow];
            tempstructT.tower_iL        = _tower_LFHCAL_iL[itow];
            tempstructT.tower_trueID    = _tower_LFHCAL_trueID[itow];
            input_towers_LFHCAL.push_back(tempstructT);
            input_towersForMC_LFHCAL.push_back(tempstructT);
            if (verbosity > 1 ) std::cout << tempstructT.tower_E << "\t" << tempstructT.tower_iEta << "\t" << tempstructT.tower_iPhi << std::endl;
        }
    }

    // sort vector in descending energy order
    std::sort(input_towers_LFHCAL.begin(), input_towers_LFHCAL.end(), &acompare);
    std::sort(input_towersForMC_LFHCAL.begin(), input_towersForMC_LFHCAL.end(), &acompareTrueID);
    
    std::vector<int> clslabels_LFHCAL;
    while (!input_towers_LFHCAL.empty()) {
        cluster_towers_LFHCAL.clear();
        clslabels_LFHCAL.clear();

        clustersStrct tempstructC;

        if(input_towers_LFHCAL.at(0).tower_E > seedEC_LFHCAL){

        if(clusterizerEnum==kC3){
            tempstructC = findCircularCluster(3, seedEC_LFHCAL, caloID, input_towers_LFHCAL, cluster_towers_LFHCAL, clslabels_LFHCAL);
        } else if (clusterizerEnum==kC5){
            tempstructC = findCircularCluster(5, seedEC_LFHCAL, caloID, input_towers_LFHCAL, cluster_towers_LFHCAL, clslabels_LFHCAL);
        } else if (clusterizerEnum==k3x3){
            tempstructC = findSquareCluster(3, seedEC_LFHCAL, caloID, input_towers_LFHCAL, cluster_towers_LFHCAL, clslabels_LFHCAL);
        } else if (clusterizerEnum==k5x5){
            tempstructC = findSquareCluster(5, seedEC_LFHCAL, caloID, input_towers_LFHCAL, cluster_towers_LFHCAL, clslabels_LFHCAL);
        } else if (clusterizerEnum==kV1){  
            tempstructC = findV1Cluster(seedEC_LFHCAL, caloID, input_towers_LFHCAL, cluster_towers_LFHCAL, clslabels_LFHCAL);
        } else if (clusterizerEnum==kV3){  
            tempstructC = findV3Cluster(seedEC_LFHCAL, caloID, input_towers_LFHCAL, cluster_towers_LFHCAL, clslabels_LFHCAL);
        } else if (clusterizerEnum==kMA){  
            tempstructC = findMACluster(seedEC_LFHCAL, caloID, input_towers_LFHCAL, cluster_towers_LFHCAL, clslabels_LFHCAL);
        } else {
            std::cout << "incorrect clusterizer selected! Enum " << clusterizerEnum << " not defined!" << std::endl;
            return;
        }
        tempstructC.cluster_towers = cluster_towers_LFHCAL;

        rec_clusters_LFHCAL.push_back(tempstructC);
        nclusters++;
      
        } else {
          input_towers_LFHCAL.clear();
        }
    }

    float firstClusterE_LFHCAL     = 0;
    if(nclusters==0) 
      firstClusterE_LFHCAL   = 0;
    else 
      firstClusterE_LFHCAL   = rec_clusters_LFHCAL.at(0).cluster_E;

    float maxSeedEnergy_LFHCAL = 0;
    for (int it = 0; it < (int)rec_clusters_LFHCAL.size(); it++){
        cout << "filling hist for cluster: "  << it << "\t E: " << rec_clusters_LFHCAL.at(it).cluster_E << "\t seed E: " << rec_clusters_LFHCAL.at(it).cluster_seed <<  endl;
        for(int tic=0; tic < (int)(rec_clusters_LFHCAL.at(it).cluster_towers).size(); tic++){
            int iEtaTwr = rec_clusters_LFHCAL.at(it).cluster_towers.at(tic).tower_iEta;
            int iPhiTwr = rec_clusters_LFHCAL.at(it).cluster_towers.at(tic).tower_iPhi;
            int iLTwr   = rec_clusters_LFHCAL.at(it).cluster_towers.at(tic).tower_iL;
            float eTwr  = rec_clusters_LFHCAL.at(it).cluster_towers.at(tic).tower_E;
            if (verbosity > 0) std::cout << "LFHCAL: " << "\t" << Form("%.3f", eTwr) << "\t eta \t" << iEtaTwr << "\t phi \t" << iPhiTwr << std::endl;
            int tmpX = arr_LFHCALGeoStrt[iEtaTwr][iPhiTwr][iLTwr].tower_x;
            int tmpY = arr_LFHCALGeoStrt[iEtaTwr][iPhiTwr][iLTwr].tower_y;
            int tmpZ = arr_LFHCALGeoStrt[iEtaTwr][iPhiTwr][iLTwr].tower_z;
                h_XYZMapEvt_Cl[it]->Fill(tmpZ+3*dZ, tmpX+dX, tmpY+dY, eTwr);
            }
        if( rec_clusters_LFHCAL.at(it).cluster_seed > maxSeedEnergy_LFHCAL )
            maxSeedEnergy_LFHCAL   = rec_clusters_LFHCAL.at(it).cluster_seed;
    }

    cout << "here" << endl;
    int nMCCurr   = 0;
    int nMCCurrID = input_towersForMC_LFHCAL.at(0).tower_trueID;
    std::vector<towersStrct> cluster_towersMC;
    std::vector<clustersStrct> MC_clusters_LFHCAL;
    clustersStrct tempstructMC;
    tempstructMC.cluster_trueID = nMCCurrID;

    while (!input_towersForMC_LFHCAL.empty() ) {
      if (nMCCurrID != input_towersForMC_LFHCAL.at(0).tower_trueID ){
            tempstructMC.cluster_Eta    = _mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_Phi    = _mcpart_Phi[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_trueID = nMCCurrID;
            tempstructMC.cluster_NTowers  = cluster_towersMC.size();
            tempstructMC.cluster_towers   = cluster_towersMC;
            MC_clusters_LFHCAL.push_back(tempstructMC);
            nMCCurrID = input_towersForMC_LFHCAL.at(0).tower_trueID;
            // reseting tempstructMC
            tempstructMC  = clustersStrct();
            tempstructMC.cluster_trueID = nMCCurrID;
            cluster_towersMC.clear();
      }
      cluster_towersMC.push_back(input_towersForMC_LFHCAL.at(0));
      tempstructMC.cluster_E += input_towersForMC_LFHCAL.at(0).tower_E;

      if (input_towersForMC_LFHCAL.size() == 1){
            tempstructMC.cluster_Eta    = _mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_Phi    = _mcpart_Phi[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_trueID = nMCCurrID;
            tempstructMC.cluster_NTowers  = cluster_towersMC.size();
            tempstructMC.cluster_towers   = cluster_towersMC;
            MC_clusters_LFHCAL.push_back(tempstructMC);    
      }
      input_towersForMC_LFHCAL.erase(input_towersForMC_LFHCAL.begin());
    }
    std::sort(MC_clusters_LFHCAL.begin(), MC_clusters_LFHCAL.end(), &acompareCl);
    cout << __LINE__ << endl;

    for (int nMC = 0; nMC < (int)MC_clusters_LFHCAL.size() && nMC < 50; nMC++  ){
        int mcID    = MC_clusters_LFHCAL.at(nMC).cluster_trueID;
        TString tempMCLabel = Form("%i: r. %.1f GeV, #eta = %1.2f", _mcpart_PDG[GetCorrectMCArrayEntry(mcID)] , MC_clusters_LFHCAL.at(nMC).cluster_E,  MC_clusters_LFHCAL.at(nMC).cluster_Eta  );
        mcPartEnergies.push_back(tempMCLabel);
        cout << tempMCLabel.Data() << endl;
        for (int tic = 0; tic < (int)(MC_clusters_LFHCAL.at(nMC).cluster_towers).size(); tic++){
            int iEtaTwr = MC_clusters_LFHCAL.at(nMC).cluster_towers.at(tic).tower_iEta;
            int iPhiTwr = MC_clusters_LFHCAL.at(nMC).cluster_towers.at(tic).tower_iPhi;
            int iLTwr   = MC_clusters_LFHCAL.at(nMC).cluster_towers.at(tic).tower_iL;
            float eTwr  = MC_clusters_LFHCAL.at(nMC).cluster_towers.at(tic).tower_E;
            int tmpX = arr_LFHCALGeoStrt[iEtaTwr][iPhiTwr][iLTwr].tower_x;
            int tmpY = arr_LFHCALGeoStrt[iEtaTwr][iPhiTwr][iLTwr].tower_y;
            int tmpZ = arr_LFHCALGeoStrt[iEtaTwr][iPhiTwr][iLTwr].tower_z;
            h_XYZMapEvt_MCPart[nMC]->Fill(tmpZ+3*dZ, tmpX+dX, tmpY+dY, eTwr);
        }  
        nMCCurr++;
    }

    // **********************************************************************************************************************************************************
    // FEMC clusterizer

    // vector of towers used in the clusterizer
    std::vector<towersStrct> input_towers_FEMC;
    std::vector<towersStrct> input_towersForMC_FEMC;
    // vector of towers within the currently found cluster
    std::vector<towersStrct> cluster_towers_FEMC;
    std::vector<clustersStrct> rec_clusters_FEMC;
    std::vector<float> recclusterEnergies_FEMC;

    caloID = 1;     // clusterization for FEMC

    // fill vector with towers for clusterization above aggregation threshold
    for(int itow=0; itow<(int)_nTowers_FEMC; itow++){
        if(_tower_FEMC_E[itow]>aggEC_FEMC){
            towersStrct tempstructT;
            tempstructT.tower_E         = _tower_FEMC_E[itow];
            tempstructT.tower_iEta      = _tower_FEMC_iEta[itow];
            tempstructT.tower_iPhi      = _tower_FEMC_iPhi[itow];
            tempstructT.tower_iL        = -1;
            tempstructT.tower_trueID    = _tower_FEMC_trueID[itow];
            input_towers_FEMC.push_back(tempstructT);
            input_towersForMC_FEMC.push_back(tempstructT);
            if (verbosity > 1 ) std::cout << tempstructT.tower_E << "\t" << tempstructT.tower_iEta << "\t" << tempstructT.tower_iPhi << std::endl;
        }
    }

    // sort vector in descending energy order
    std::sort(input_towers_FEMC.begin(), input_towers_FEMC.end(), &acompare);
    std::sort(input_towersForMC_FEMC.begin(), input_towersForMC_FEMC.end(), &acompareTrueID);
    
    std::vector<int> clslabels_FEMC;
    while (!input_towers_FEMC.empty()) {
        cluster_towers_FEMC.clear();
        clslabels_FEMC.clear();

        clustersStrct tempstructC;

        if(input_towers_FEMC.at(0).tower_E > seedEC_FEMC){

        if(clusterizerEnum==kC3){
            tempstructC = findCircularCluster(3, seedEC_FEMC, caloID, input_towers_FEMC, cluster_towers_FEMC, clslabels_FEMC);
        } else if (clusterizerEnum==kC5){
            tempstructC = findCircularCluster(5, seedEC_FEMC, caloID, input_towers_FEMC, cluster_towers_FEMC, clslabels_FEMC);
        } else if (clusterizerEnum==k3x3){
            tempstructC = findSquareCluster(3, seedEC_FEMC, caloID, input_towers_FEMC, cluster_towers_FEMC, clslabels_FEMC);
        } else if (clusterizerEnum==k5x5){
            tempstructC = findSquareCluster(5, seedEC_FEMC, caloID, input_towers_FEMC, cluster_towers_FEMC, clslabels_FEMC);
        } else if (clusterizerEnum==kV1){  
            tempstructC = findV1Cluster(seedEC_FEMC, caloID, input_towers_FEMC, cluster_towers_FEMC, clslabels_FEMC);
        } else if (clusterizerEnum==kV3){  
            tempstructC = findV3Cluster(seedEC_FEMC, caloID, input_towers_FEMC, cluster_towers_FEMC, clslabels_FEMC);
        } else if (clusterizerEnum==kMA){  
            tempstructC = findMACluster(seedEC_FEMC, caloID, input_towers_FEMC, cluster_towers_FEMC, clslabels_FEMC);
        } else {
            std::cout << "incorrect clusterizer selected! Enum " << clusterizerEnum << " not defined!" << std::endl;
            return;
        }
        tempstructC.cluster_towers = cluster_towers_FEMC;
        
        rec_clusters_FEMC.push_back(tempstructC);
        nclusters++;
      
        } else {
          input_towers_FEMC.clear();
        }
    }

    float firstClusterE_FEMC     = 0;
    if(nclusters==0) 
      firstClusterE_FEMC   = 0;
    else 
      firstClusterE_FEMC   = rec_clusters_FEMC.at(0).cluster_E;

    float maxSeedEnergy_FEMC = 0;
    for( int it=0; it < (int)rec_clusters_FEMC.size(); it++){
        cout << "filling hist for cluster: "  << (int)rec_clusters_LFHCAL.size()+it << "\t E: " << rec_clusters_FEMC.at(it).cluster_E << "\t seed E: " << rec_clusters_FEMC.at(it).cluster_seed <<  endl;
        cout <<(int)rec_clusters_FEMC.size() << endl;
        for (int tic = 0; tic < (int)rec_clusters_FEMC.at(it).cluster_towers.size(); tic++){
            int iEtaTwr = rec_clusters_FEMC.at(it).cluster_towers.at(tic).tower_iEta;
            int iPhiTwr = rec_clusters_FEMC.at(it).cluster_towers.at(tic).tower_iPhi;
            float eTwr  = rec_clusters_FEMC.at(it).cluster_towers.at(tic).tower_E;
            if (verbosity > 0) std::cout << "FEMC:" << "\t" << Form("%.3f",eTwr) << "\t eta \t" << iEtaTwr << "\t phi \t" << iPhiTwr << std::endl;
            int tmpX = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_x;
            int tmpY = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_y;
            int tmpZ = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_z;
            h_XYZMapEvt_Cl[(int)rec_clusters_LFHCAL.size()+it]->Fill(tmpZ+7*dZ, tmpX, tmpY, eTwr);
        }
    }

    nMCCurrID = input_towersForMC_FEMC.at(0).tower_trueID;
    cluster_towersMC.clear();
    std::vector<clustersStrct> MC_clusters;
    tempstructMC.cluster_trueID = nMCCurrID;

    while (!input_towersForMC_FEMC.empty() ) {
      if (nMCCurrID != input_towersForMC_FEMC.at(0).tower_trueID ){
            tempstructMC.cluster_Eta    = _mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_Phi    = _mcpart_Phi[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_trueID = nMCCurrID;
            tempstructMC.cluster_NTowers  = cluster_towersMC.size();
            tempstructMC.cluster_towers   = cluster_towersMC;
            MC_clusters.push_back(tempstructMC);
            nMCCurrID = input_towersForMC_FEMC.at(0).tower_trueID;
            // reseting tempstructMC
            tempstructMC  = clustersStrct();
            tempstructMC.cluster_trueID = nMCCurrID;
            cluster_towersMC.clear();
      }
      cluster_towersMC.push_back(input_towersForMC_FEMC.at(0));
      tempstructMC.cluster_E += input_towersForMC_FEMC.at(0).tower_E;

      if (input_towersForMC_LFHCAL.size() == 1){
            tempstructMC.cluster_Eta    = _mcpart_Eta[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_Phi    = _mcpart_Phi[GetCorrectMCArrayEntry(nMCCurrID)] ;
            tempstructMC.cluster_trueID = nMCCurrID;
            tempstructMC.cluster_NTowers  = cluster_towersMC.size();
            tempstructMC.cluster_towers   = cluster_towersMC;
            MC_clusters.push_back(tempstructMC);    
      }
      input_towersForMC_FEMC.erase(input_towersForMC_FEMC.begin());
    }
    std::sort(MC_clusters.begin(), MC_clusters.end(), &acompareCl);
    cout << __LINE__ << endl;


    std::vector<Bool_t> checkedClusters_FEMC((int)MC_clusters.size(), kFALSE);
    for( int nLF = 0; nLF < nMCCurr; nLF++){
        int tmpTrueId = MC_clusters_LFHCAL.at(nLF).cluster_trueID;
        for(int nMC = 0; nMC < (int) MC_clusters.size() && nMC < 50; nMC++){
            int mcID = MC_clusters.at(nMC).cluster_trueID;
            if( mcID == tmpTrueId ){
                for (int tic = 0; tic < (int)(MC_clusters.at(nMC).cluster_towers).size(); tic++){
                    int iEtaTwr = MC_clusters.at(nMC).cluster_towers.at(tic).tower_iEta;
                    int iPhiTwr = MC_clusters.at(nMC).cluster_towers.at(tic).tower_iPhi;
                    float eTwr  = MC_clusters.at(nMC).cluster_towers.at(tic).tower_E;
                    int tmpX = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_x;
                    int tmpY = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_y;
                    int tmpZ = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_z;
                    h_XYZMapEvt_MCPart[nLF]->Fill(tmpZ+7*dZ, tmpX, tmpY, eTwr);
                }
            checkedClusters_FEMC.at(nMC) = kTRUE;
            }
        }
    }
    
    for (int nMC = 0; nMC < (int)MC_clusters.size() && nMC < 50; nMC++  ){
        if( checkedClusters_FEMC.at(nMC) ) continue;
        int mcID    = MC_clusters.at(nMC).cluster_trueID;
        TString tempMCLabel = Form("%i: r. %.1f GeV, #eta = %1.2f", _mcpart_PDG[GetCorrectMCArrayEntry(mcID)] , MC_clusters.at(nMC).cluster_E,  MC_clusters.at(nMC).cluster_Eta  );
        mcPartEnergies.push_back(tempMCLabel);
        for (int tic = 0; tic < (int)(MC_clusters.at(nMC).cluster_towers).size(); tic++){
            int iEtaTwr = MC_clusters.at(nMC).cluster_towers.at(tic).tower_iEta;
            int iPhiTwr = MC_clusters.at(nMC).cluster_towers.at(tic).tower_iPhi;
            float eTwr  = MC_clusters.at(nMC).cluster_towers.at(tic).tower_E;
            int tmpX = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_x;
            int tmpY = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_y;
            int tmpZ = arr_FEMCGeoStrt[iEtaTwr][iPhiTwr].tower_z;
            h_XYZMapEvt_MCPart[nMCCurr]->Fill(tmpZ+7*dZ, tmpX, tmpY, eTwr);
        }  
        nMCCurr++;
    }

    if (sumedE_LFHCAL < 0.01 && sumedE_FEMC < 0.01){
        cout << "total E calo LFHCAL: " << sumedE_LFHCAL << " total E calo FEMC: " << sumedE_FEMC << " No particle hit the calorimeter! Aborting." << endl; 
        return;
    } else {
        cout << "total E calo LFHCAL: " << sumedE_LFHCAL << " total E calo FEMC: " << sumedE_FEMC << endl; 
    }

    float maxSeedEnergy = maxSeedEnergy_LFHCAL;
    if( maxSeedEnergy_FEMC > maxSeedEnergy ) maxSeedEnergy = maxSeedEnergy_FEMC;

    // ********************************************************************************************************************************************************
    // PLOTTING

    TCanvas* cReso = new TCanvas("cReso","",0,0,1000,950);
    DrawGammaCanvasSettings( cReso, 0.12, 0.05, 0.01, 0.1);
    cReso->SetTheta(16);
    cReso->SetPhi(68);

    labelDet += " FEMC";

    // *****************************************************************************************
    // 3D hist with box option - all hits
    //******************************************************************************************

    SetStyleHistoTH3ForGraphs(h_XYZMapEvt_Hole, "tower z index", "tower x index","tower y index", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1, 1.1, 1.15, 505, 510,510);
    h_XYZMapEvt_Hole->SetLineColor(kGray);
    h_XYZMapEvt_Hole->Draw("box");
    h_XYZMapEvt_Ecut->SetMaximum(maxSeedEnergy);
    h_XYZMapEvt_Ecut->SetLineColor(kGray+1);
    h_XYZMapEvt_Ecut->Draw("same,box");

    drawLatexAdd(labelEnergy.Data(),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(labelDet,0.02,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print( Form("%s/3D/3D_Event%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()) );
    std::cout <<  Form("%s/3D/3D_Event%d_Ecut.%s", outputDir.Data(),plotEvent, suffix.Data()) << std::endl;

    // *****************************************************************************************
    // 3D hist with box option - clusters
    //******************************************************************************************
    cout << "Nclusters: " << nclusters << endl;
    for(int nCl=0; nCl < nclusters; nCl++){
        h_XYZMapEvt_Cl[nCl]->SetMaximum(maxSeedEnergy);
        h_XYZMapEvt_Cl[nCl]->SetLineColor(colorCluster[nCl]);
        h_XYZMapEvt_Cl[nCl]->Draw("same,box");
        cout << nCl << endl;
    }

    drawLatexAdd(labelEnergy.Data(),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(labelDet,0.02,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    std::cout << Form("%s/3DwClusters/3D_Event%d_Ecut_clusters.%s", outputDir.Data(),plotEvent, suffix.Data()) << std::endl;
    cReso->Print(Form("%s/3DwClusters/3D_Event%d_Ecut_clusters.%s", outputDir.Data(),plotEvent, suffix.Data()));


    // *****************************************************************************************
    // 3D hist with box option and different MC particles ontop
    //******************************************************************************************
    h_XYZMapEvt_Hole->Draw("box");
    h_XYZMapEvt_Ecut->SetMaximum(maxSeedEnergy);
    h_XYZMapEvt_Ecut->SetLineColor(kGray+2);
    h_XYZMapEvt_Ecut->Draw("same,box");

    cout << "nMCCurr: " << nMCCurr << endl;
    for(int nMC=0; nMC < nMCCurr; nMC++){
        if( h_XYZMapEvt_MCPart[nMC] ){
            h_XYZMapEvt_MCPart[nMC]->SetLineColor(colorCluster[nMC]);
            h_XYZMapEvt_MCPart[nMC]->Draw("box,same");
            int cols = nMC%2;
            int rows = (nMC/2);
            drawLatexAdd((TString)(mcPartEnergies.at(nMC)).Data() , 0.44+0.27*cols, 0.93-(rows*0.55*textSizeLabelsRel), 0.55*textSizeLabelsRel, kFALSE, kFALSE, kFALSE, colorCluster[nMC]);
        }
    }

    drawLatexAdd(labelEnergy.Data(),0.02,0.95,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(labelDet,0.02,0.90,textSizeLabelsRel,kFALSE,kFALSE,kFALSE);

    cReso->Print( Form("%s/3DwClusters/3D_Event%d_Ecut_MCPart.%s", outputDir.Data(),plotEvent, suffix.Data()) );
    std::cout <<  Form("%s/3DwClusters/3D_Event%d_Ecut_MCPart.%s", outputDir.Data(),plotEvent, suffix.Data()) << std::endl;

}

