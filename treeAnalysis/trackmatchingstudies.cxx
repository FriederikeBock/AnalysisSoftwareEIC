#include <algorithm>
// ANCHOR debug output verbosity
int verbosityTM = 1;


TH2F*  h_TMstudies_x_y_clus[_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_x_y_true[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_x_y_track[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_x_y_projection[_active_calo]; // [calorimeter_enum][algorithm_enum]

TH1F*  h_TMstudies_clusSpec[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_Eta[_active_calo][nEta+1]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_Matched[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_Matched_Eta[_active_calo][nEta+1]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_MatchedCalo[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_MatchedCalo_EDiff1[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_MatchedCalo_EDiff2[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_special[_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_clusSpec_special_Matched[_active_calo]; // [calorimeter_enum][algorithm_enum]

const int diffCases_TMStudies = 1;
TString strCases_TMStudies[diffCases_TMStudies] = {"all"};//, "leading"};

TH2F*  h_TMstudies_dx_dy_track[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_dy_projection[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_dy_TTLhits[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_dy_calomatch[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_E1_E2_calomatch[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsClsEta_track[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsClsEta_projection[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsTrkEta_track[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_dx_vsTrkEta_projection[diffCases_TMStudies][_active_calo][_active_algo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_isMatched[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_E[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_PDG[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH1F*  h_TMstudies_isMatched_findable[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_findable_E[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
TH2F*  h_TMstudies_isMatched_findable_PDG[diffCases_TMStudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
// const in iPDG_TMstudies =4;
// int iTMstudies_partPDG[iPDG_TMstudies] = {11,211,321,2212};
// TString strTMstudies_partPDG[iPDG_TMstudies] = {"electron","cpion","ckaon","proton"};
// TH2F*  h_TMstudies_isMatchedCell_track[iPDG_TMstudies][_active_calo]; // [calorimeter_enum][algorithm_enum]
// TH2F*  h_TMstudies_isMatchedCell_projection[iPDG_TMstudies][_active_calo]; // [calorimeter_enum][algorithm_enum]

  int _nclusters_TM = 0;
  float* clusters_TM_E            = new float[_maxNclusters];
  float* clusters_TM_Eta         = new float[_maxNclusters];
  float* clusters_TM_Phi         = new float[_maxNclusters];
  float* clusters_TM_M02         = new float[_maxNclusters];
  float* clusters_TM_M20         = new float[_maxNclusters];
  bool* clusters_TM_isMatched         = new bool[_maxNclusters];
  int* clusters_TM_NTower         = new int[_maxNclusters];
  int* clusters_TM_trueID       = new int[_maxNclusters];
  int* clusters_TM_NtrueID       = new int[_maxNclusters];
  float* clusters_TM_X         = new float[_maxNclusters];
  float* clusters_TM_Y         = new float[_maxNclusters];
  float* clusters_TM_Z         = new float[_maxNclusters];

    int nelecinaccFHCAL_TMstudies = 0;
    int npioninaccFHCAL_TMstudies = 0;
    int partinaccFHCAL_TMstudies = 0;

    int nelecinaccFEMC_TMstudies = 0;
    int npioninaccFEMC_TMstudies = 0;
    int partinaccFEMC_TMstudies = 0;

// ANCHOR track/projection matching function
// TODO at some point we might want E/p
bool trackmatchingstudies(){
  bool b_TMstudies_x_y_track_filled[_active_calo] = {false};
  for(int icalo=0;icalo<_active_calo;icalo++){
      if(!h_TMstudies_x_y_projection[icalo])h_TMstudies_x_y_projection[icalo] 	= new TH2F(Form("h_TMstudies_x_y_projection_%s",str_calorimeter[icalo].Data()), "", 299,-299,299, 299,-299,299);
      if(!h_TMstudies_x_y_track[icalo])h_TMstudies_x_y_track[icalo] 	= new TH2F(Form("h_TMstudies_x_y_track_%s",str_calorimeter[icalo].Data()), "", 299,-299,299, 299,-299,299);
      if(!h_TMstudies_x_y_true[icalo])h_TMstudies_x_y_true[icalo] 	= new TH2F(Form("h_TMstudies_x_y_true_%s",str_calorimeter[icalo].Data()), "", 299,-299,299, 299,-299,299);
      if(!h_TMstudies_clusSpec[icalo])h_TMstudies_clusSpec[icalo] 	= new TH1F(Form("h_TMstudies_clusSpec_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      if(!h_TMstudies_clusSpec_Matched[icalo])h_TMstudies_clusSpec_Matched[icalo] 	= new TH1F(Form("h_TMstudies_clusSpec_Matched_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      if(!h_TMstudies_clusSpec_MatchedCalo[icalo])h_TMstudies_clusSpec_MatchedCalo[icalo] 	= new TH1F(Form("h_TMstudies_clusSpec_MatchedCalo_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      if(!h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo])h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo] 	= new TH1F(Form("h_TMstudies_clusSpec_MatchedCalo_EDiff1_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      if(!h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo])h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo] 	= new TH1F(Form("h_TMstudies_clusSpec_MatchedCalo_EDiff2_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      if(!h_TMstudies_clusSpec_special[icalo])h_TMstudies_clusSpec_special[icalo] 	= new TH1F(Form("h_TMstudies_clusSpec_special_%s",str_calorimeter[icalo].Data()), "", 120,0,60);
      if(!h_TMstudies_clusSpec_special_Matched[icalo])h_TMstudies_clusSpec_special_Matched[icalo] 	= new TH1F(Form("h_TMstudies_clusSpec_special_Matched_%s",str_calorimeter[icalo].Data()), "", 120,0,60);

      for (Int_t et = 0; et<nEta+1; et++){
        Double_t etaMin = partEta[0];
        Double_t etaMax = partEta[nEta];
        if (et < nEta){
          etaMin = partEta[et];
          etaMax = partEta[et+1];
        }
        if(!h_TMstudies_clusSpec_Eta[icalo][et])h_TMstudies_clusSpec_Eta[icalo][et] 	= new TH1F(Form("h_TMstudies_clusSpec_Eta_%1.1f_%1.1f_%s",etaMin,etaMax,str_calorimeter[icalo].Data()), "", 120,0,60);
        if(!h_TMstudies_clusSpec_Matched_Eta[icalo][et])h_TMstudies_clusSpec_Matched_Eta[icalo][et] 	= new TH1F(Form("h_TMstudies_clusSpec_Matched_Eta_%1.1f_%1.1f_%s",etaMin,etaMax,str_calorimeter[icalo].Data()), "", 120,0,60);
      }
  }
  for(int icalo=0;icalo<_active_calo;icalo++){
    bool b_TMstudies_x_y_projection_filled = false;

    for(int icase=0;icase<diffCases_TMStudies;icase++){
      if(!h_TMstudies_isMatched[icase][icalo]){
        h_TMstudies_isMatched[icase][icalo] 	= new TH1F(Form("h_TMstudies_%s_isMatched_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*4,0,_active_algo*4);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
          h_TMstudies_isMatched[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+4, "hit matched");
        }
      }
      if(!h_TMstudies_isMatched_PDG[icase][icalo]){
        h_TMstudies_isMatched_PDG[icase][icalo] 	= new TH2F(Form("h_TMstudies_%s_isMatched_PDG_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*4,0,_active_algo*4,2500,0,2500);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
          h_TMstudies_isMatched_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+4, "hit matched");
        }
      }
      if(!h_TMstudies_isMatched_E[icase][icalo]){
        h_TMstudies_isMatched_E[icase][icalo] 	= new TH2F(Form("h_TMstudies_%s_isMatched_E_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*4,0,_active_algo*4,100,0,100);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
          h_TMstudies_isMatched_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+4, "hit matched");
        }
      }

      if(!h_TMstudies_isMatched_findable[icase][icalo]){
        h_TMstudies_isMatched_findable[icase][icalo] 	= new TH1F(Form("h_TMstudies_%s_isMatched_findable_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*4,0,_active_algo*4);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
          h_TMstudies_isMatched_findable[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+4, "hit matched");
        }
      }
      if(!h_TMstudies_isMatched_findable_E[icase][icalo]){
        h_TMstudies_isMatched_findable_E[icase][icalo] 	= new TH2F(Form("h_TMstudies_%s_isMatched_findable_E_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*4,0,_active_algo*4,100,0,100);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
          h_TMstudies_isMatched_findable_E[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+4, "hit matched");
        }
      }
      if(!h_TMstudies_isMatched_findable_PDG[icase][icalo]){
        h_TMstudies_isMatched_findable_PDG[icase][icalo] 	= new TH2F(Form("h_TMstudies_%s_isMatched_findable_PDG_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", _active_algo*4,0,_active_algo*4,2500,0,2500);
        for(int ialgo=0;ialgo<_active_algo;ialgo++){
          h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+1, Form("%s clusters", str_clusterizer[ialgo].Data()));
          h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+2, "track matched");
          h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+3, "proj matched");
          h_TMstudies_isMatched_findable_PDG[icase][icalo]->GetXaxis()->SetBinLabel(ialgo*3+4, "hit matched");
        }
      }
    }
    for(int imc=0; imc<_nMCPart; imc++){
      TVector3 mcpartvec(_mcpart_px[imc],_mcpart_py[imc],_mcpart_pz[imc]);
      float zHC = 400;
      if(icalo==kFHCAL)zHC=340;
      if(icalo==kFEMC)zHC=287;
      mcpartvec*=(zHC/mcpartvec.Z());
      // float trueeta = mcpartvec.Eta();
      // if(trueeta<3.5 && trueeta > 1.1){
        if(icalo==kFHCAL && TMath::Sqrt(TMath::Power(mcpartvec.X(),2)+TMath::Power(mcpartvec.Y(),2))<262 && TMath::Sqrt(TMath::Power(mcpartvec.Y(),2)+TMath::Power(mcpartvec.X()-6,2))>30.){
          h_TMstudies_x_y_true[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
          if(abs(_mcpart_PDG[imc])==11)nelecinaccFHCAL_TMstudies++;
          if(abs(_mcpart_PDG[imc])==211)npioninaccFHCAL_TMstudies++;
          partinaccFHCAL_TMstudies++;
        }
        if(icalo==kFEMC && TMath::Sqrt(TMath::Power(mcpartvec.X(),2)+TMath::Power(mcpartvec.Y(),2))<183 && TMath::Sqrt(TMath::Power(mcpartvec.Y(),2)+TMath::Power(mcpartvec.X()-6,2))>22.14){
          h_TMstudies_x_y_true[icalo]->Fill(mcpartvec.X(),mcpartvec.Y());
          if(abs(_mcpart_PDG[imc])==11)nelecinaccFEMC_TMstudies++;
          if(abs(_mcpart_PDG[imc])==211)npioninaccFEMC_TMstudies++;
          partinaccFEMC_TMstudies++;
        }
        // }
      // float truept  = mcpartvec.Pt();
      // float truep   = mcpartvec.Mag();

    }
    for(int ialgo=0;ialgo<_active_algo;ialgo++){
      // for(int iPDG=0;iPDG<iPDG_TMstudies;iPDG++){
      //   if(!h_TMstudies_isMatchedCell_track[iPDG][icalo])h_TMstudies_isMatchedCell_track[iPDG][icalo] 	= new TH2F(Form("h_TMstudies_isMatchedCell_track_%s_%s",strTMstudies_partPDG[iPDG].Data(),strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", 78,-39,39, 78,-39,39);
      //   if(!h_TMstudies_isMatchedCell_projection[iPDG][icalo])h_TMstudies_isMatchedCell_projection[iPDG][icalo] 	= new TH2F(Form("h_TMstudies_isMatchedCell_projection_%s_%s",strTMstudies_partPDG[iPDG].Data(),strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data()), "", 78,-39,39, 78,-39,39);
      // }
      for(int icase=0;icase<diffCases_TMStudies;icase++){
        if(!h_TMstudies_dx_dy_track[icase][icalo][ialgo])h_TMstudies_dx_dy_track[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_dy_track_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 78,-39,39);
        if(!h_TMstudies_dx_dy_projection[icase][icalo][ialgo])h_TMstudies_dx_dy_projection[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_dy_projection_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 78,-39,39);
        if(!h_TMstudies_dx_dy_TTLhits[icase][icalo][ialgo])h_TMstudies_dx_dy_TTLhits[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_dy_TTLhits_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 78,-39,39);
        if(!h_TMstudies_dx_dy_calomatch[icase][icalo][ialgo])h_TMstudies_dx_dy_calomatch[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_dy_calomatch_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 78,-39,39);
        if(!h_TMstudies_E1_E2_calomatch[icase][icalo][ialgo])h_TMstudies_E1_E2_calomatch[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_E1_E2_calomatch_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 100,0, 50.,100,0, 50.);
        if(!h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo])h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_vsClsEta_track_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 40,1.,4.);
        if(!h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo])h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_vsClsEta_projection_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 40,1.,4.);
        if(!h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo])h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_vsTrkEta_track_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 40,1.,4.);
        if(!h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo])h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo] 	= new TH2F(Form("h_TMstudies_%s_dx_vsTrkEta_projection_%s_%s",strCases_TMStudies[icase].Data(),str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 78,-39,39, 40,1.,4.);
      }
      if(!h_TMstudies_x_y_clus[icalo][ialgo])h_TMstudies_x_y_clus[icalo][ialgo] 	= new TH2F(Form("h_TMstudies_x_y_clus_%s_%s",str_calorimeter[icalo].Data(),str_clusterizer[ialgo].Data()), "", 299,-299,299, 299,-299,299);

      if(!loadClusterizerInput(
        ialgo,
        icalo,
        _nclusters_TM,
        clusters_TM_E,
        clusters_TM_Eta,
        clusters_TM_Phi,
        clusters_TM_M02,
        clusters_TM_M20,
        clusters_TM_isMatched,
        clusters_TM_NTower,
        clusters_TM_trueID,
        clusters_TM_NtrueID,
        clusters_TM_X,
        clusters_TM_Y,
        clusters_TM_Z))continue;

      if(_nclusters_TM==0) continue;
      int projectionlayer = ReturnProjectionIndexForCalorimeter(icalo);

      float zHC = 400;
      if(icalo==kFHCAL)zHC=340;
      if(icalo==kFEMC)zHC=287;
      // TODO MUST OPTIMIZE HERE, MAYBE LOAD TOWER POSITION FROM GEOMETRY TREE????
      // TODO MUST OPTIMIZE HERE, MAYBE LOAD TOWER POSITION FROM GEOMETRY TREE????
      // TODO MUST OPTIMIZE HERE, MAYBE LOAD TOWER POSITION FROM GEOMETRY TREE????
      // TODO MUST OPTIMIZE HERE, MAYBE LOAD TOWER POSITION FROM GEOMETRY TREE????

      float matchingwdw = 10.;
      if(icalo==kFEMC) matchingwdw = 5.5;
      // float matchingwdw = 15.;
      // if(icalo==kFEMC) matchingwdw = 8.;

      std::vector<int> usedIDs;
      std::vector<int> currentIDs;
      std::vector<int> usedMClabels;
      for(Int_t iclus=0; iclus<_nclusters_TM; iclus++){
        currentIDs.clear();
        bool ismatched_track = false;
        bool ismatched_projection = false;
        bool ismatched_hit = false;
        bool ismatched_largest_track = false;
        bool ismatched_largest_projection = false;
        if((clusters_TM_M02[iclus]<0.1))continue;
        h_TMstudies_clusSpec[icalo]->Fill(clusters_TM_E[iclus]);
        if((icalo==kFHCAL && abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==211) || (icalo==kFEMC && abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==11)){
          h_TMstudies_clusSpec_special[icalo]->Fill(clusters_TM_E[iclus]);
          for (Int_t et = 0; et<nEta+1; et++){
            Double_t etaMin = partEta[0];
            Double_t etaMax = partEta[nEta];
            if (et < nEta){
              etaMin = partEta[et];
              etaMax = partEta[et+1];
            }
            if(clusters_TM_Eta[iclus]>etaMin && clusters_TM_Eta[iclus]<etaMax){
              h_TMstudies_clusSpec_Eta[icalo][et]->Fill(clusters_TM_E[iclus]);
            }
          }
        }
        // if(!(clusters_TM_Eta[iclus]>3.0))continue;
        // if((clusters_TM_NTower[iclus]<2))continue;
        // if((clusters_TM_Eta[iclus]>3.2))continue;

        // find largest cluster for current MC label and make matching plots for only the largest cluster
        usedIDs.push_back(iclus);
        currentIDs.push_back(iclus);
        for(Int_t icluslbl=0; icluslbl<_nclusters_TM; icluslbl++){
          if((int)(clusters_TM_trueID[iclus])==(int)(clusters_TM_trueID[icluslbl])){

            if ( !(std::find(usedIDs.begin(), usedIDs.end(), icluslbl) != usedIDs.end()) ){
              usedIDs.push_back(icluslbl);
              currentIDs.push_back(icluslbl);
            }
          }
        }

        int largestClusterID = -1;
        float largestClusterE = -1;
        for (int id = 0; id < currentIDs.size() ; id++){
          if ( clusters_TM_E[currentIDs.at(id)] > largestClusterE )
            largestClusterE = clusters_TM_E[currentIDs.at(id)];
            largestClusterID = currentIDs.at(id);
        }
        bool fillcurrentlargest = false;
        if ( diffCases_TMStudies>1 && !(std::find(usedMClabels.begin(), usedMClabels.end(), (int)(clusters_TM_trueID[iclus])) != usedMClabels.end()) ){
          usedMClabels.push_back((int)(clusters_TM_trueID[iclus]));
          fillcurrentlargest = true;
        }

        // bool matchfound = false;
        h_TMstudies_isMatched[0][icalo]->Fill(3*ialgo);
        h_TMstudies_isMatched_E[0][icalo]->Fill(3*ialgo,clusters_TM_E[iclus] );
        h_TMstudies_isMatched_PDG[0][icalo]->Fill(3*ialgo,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
        if(fillcurrentlargest){
          h_TMstudies_isMatched[1][icalo]->Fill(3*ialgo);
          h_TMstudies_isMatched_E[1][icalo]->Fill(3*ialgo,clusters_TM_E[largestClusterID] );
          h_TMstudies_isMatched_PDG[1][icalo]->Fill(3*ialgo,abs(_mcpart_PDG[clusters_TM_trueID[largestClusterID]]) );
        }
        bool isfindable = true;
        if(abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==22 || abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==2112 || abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==130 || abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==13) isfindable = false;
        if(isfindable){
          h_TMstudies_isMatched_findable[0][icalo]->Fill(3*ialgo);
          h_TMstudies_isMatched_findable_E[0][icalo]->Fill(3*ialgo,clusters_TM_E[iclus] );
          h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(3*ialgo,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
          if(fillcurrentlargest){
            h_TMstudies_isMatched_findable[1][icalo]->Fill(3*ialgo);
            h_TMstudies_isMatched_findable_E[1][icalo]->Fill(3*ialgo,clusters_TM_E[largestClusterID] );
            h_TMstudies_isMatched_findable_PDG[1][icalo]->Fill(3*ialgo,abs(_mcpart_PDG[clusters_TM_trueID[largestClusterID]]) );
          }
        }
        h_TMstudies_x_y_clus[icalo][ialgo]->Fill(clusters_TM_X[iclus],clusters_TM_Y[iclus]);
        for(Int_t iproj=0; iproj<_nProjections; iproj++){
          // if(matchfound) continue;
          // if((abs(_track_Proj_x[iclus])<15 && abs(_track_Proj_y[iclus])<15))continue;
          if(_track_Proj_t[iproj]<-9000) continue;
          if(_track_Proj_x[iproj]==0 && _track_Proj_y[iproj]==0) continue;
          if(!b_TMstudies_x_y_projection_filled){
            if(_track_ProjLayer[iproj]==2)h_TMstudies_x_y_projection[kFHCAL]->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
            if(_track_ProjLayer[iproj]==1)h_TMstudies_x_y_projection[kFEMC]->Fill(_track_Proj_x[iproj],_track_Proj_y[iproj]);
          }
          if(_track_ProjLayer[iproj]!=projectionlayer) continue;
          if(clusters_TM_X[iclus]==0 && clusters_TM_Y[iclus]==0) continue;
          // check eta difference
          TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
          float projeta = projvec.Eta();
          h_TMstudies_dx_dy_projection[0][icalo][ialgo]->Fill(_track_Proj_x[iproj]-clusters_TM_X[iclus],_track_Proj_y[iproj]-clusters_TM_Y[iclus]);
          h_TMstudies_dx_vsClsEta_projection[0][icalo][ialgo]->Fill(_track_Proj_x[iproj]-clusters_TM_X[iclus],clusters_TM_Eta[iclus]);
          h_TMstudies_dx_vsTrkEta_projection[0][icalo][ialgo]->Fill(_track_Proj_x[iproj]-clusters_TM_X[iclus],projeta);
          if(fillcurrentlargest){
            h_TMstudies_dx_dy_projection[1][icalo][ialgo]->Fill(_track_Proj_x[iproj]-clusters_TM_X[largestClusterID],_track_Proj_y[iproj]-clusters_TM_Y[largestClusterID]);
            h_TMstudies_dx_vsClsEta_projection[1][icalo][ialgo]->Fill(_track_Proj_x[iproj]-clusters_TM_X[largestClusterID],clusters_TM_Eta[largestClusterID]);
            h_TMstudies_dx_vsTrkEta_projection[1][icalo][ialgo]->Fill(_track_Proj_x[iproj]-clusters_TM_X[largestClusterID],projeta);
          }
          // if(TMath::Sqrt(TMath::Power(_track_Proj_x[iproj]-clusters_TM_X[iclus],2)+TMath::Power(_track_Proj_y[iproj]-clusters_TM_Y[iclus],2)) <15)matchfound=true;
          // TVector3 projvec(_track_Proj_x[iproj],_track_Proj_y[iproj],_track_Proj_z[iproj]);
          // float projeta = projvec.Eta();
          // float projphi = (projvec.Phi()<0 ? projvec.Phi()+TMath::Pi() : projvec.Phi()-TMath::Pi());
          // h_TMstudies_dx_dy_projection[0][icalo][ialgo]->Fill(projeta-clusters_TM_Eta[iclus],projphi-clusters_TM_Phi[iclus]);
          if(abs(_track_Proj_x[iproj]-clusters_TM_X[iclus]) < matchingwdw){
            if(abs(_track_Proj_y[iproj]-clusters_TM_Y[iclus]) < matchingwdw){
              ismatched_projection = true;
              h_TMstudies_clusSpec_Matched[icalo]->Fill(clusters_TM_E[iclus]);
              if((icalo==kFHCAL && abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==211) || (icalo==kFEMC && abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==11)){
                h_TMstudies_clusSpec_special_Matched[icalo]->Fill(clusters_TM_E[iclus]);
                for (Int_t et = 0; et<nEta+1; et++){
                  Double_t etaMin = partEta[0];
                  Double_t etaMax = partEta[nEta];
                  if (et < nEta){
                    etaMin = partEta[et];
                    etaMax = partEta[et+1];
                  }
                  if(clusters_TM_Eta[iclus]>etaMin && clusters_TM_Eta[iclus]<etaMax){
                    h_TMstudies_clusSpec_Matched_Eta[icalo][et]->Fill(clusters_TM_E[iclus]);
                  }
                }
              }
            }
          }
          if(abs(_track_Proj_x[iproj]-clusters_TM_X[largestClusterID]) < matchingwdw){
            if(abs(_track_Proj_y[iproj]-clusters_TM_Y[largestClusterID]) < matchingwdw){
              ismatched_largest_projection = true;
            }
          }
        }
        b_TMstudies_x_y_projection_filled = true;

        for(Int_t itrk=0; itrk<_nTracks; itrk++){
          // check eta difference
          TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
          trackvec*=(zHC/trackvec.Z());
          float tracketa = trackvec.Eta();
          if(!b_TMstudies_x_y_track_filled[icalo])h_TMstudies_x_y_track[icalo]->Fill(trackvec.X(),trackvec.Y());

          if(clusters_TM_X[iclus]==0 && clusters_TM_Y[iclus]==0) continue;

          h_TMstudies_dx_dy_track[0][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[iclus],trackvec.Y()-clusters_TM_Y[iclus]);
          h_TMstudies_dx_vsClsEta_track[0][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[iclus],clusters_TM_Eta[iclus]);
          h_TMstudies_dx_vsTrkEta_track[0][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[iclus],tracketa);
          if(fillcurrentlargest){
            h_TMstudies_dx_dy_track[1][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[largestClusterID],trackvec.Y()-clusters_TM_Y[largestClusterID]);
            h_TMstudies_dx_vsClsEta_track[1][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[largestClusterID],clusters_TM_Eta[largestClusterID]);
            h_TMstudies_dx_vsTrkEta_track[1][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[largestClusterID],tracketa);
          }
          if(abs(trackvec.X()-clusters_TM_X[iclus]) < matchingwdw){
            if(abs(trackvec.Y()-clusters_TM_Y[iclus]) < matchingwdw){
              ismatched_track = true;
            }
          }
          if(abs(trackvec.X()-clusters_TM_X[largestClusterID]) < matchingwdw){
            if(abs(trackvec.Y()-clusters_TM_Y[largestClusterID]) < matchingwdw){
              ismatched_largest_track = true;
            }
          }
        }
        b_TMstudies_x_y_track_filled[icalo] = true;

        for(Int_t ihit=0; ihit<_nHitsLayers; ihit++){
          if (_hits_layerID[ihit] == -1) continue;
          if( (icalo==kFHCAL && _hits_layerID[ihit]==2) ||
              (icalo==kFEMC && _hits_layerID[ihit]==0)){
            // h_hits_layer_xy[icalo]->Fill(_hits_x[ihit],_hits_y[ihit]);

            if(clusters_TM_X[iclus]==0 && clusters_TM_Y[iclus]==0) continue;
            h_TMstudies_dx_dy_TTLhits[0][icalo][ialgo]->Fill(_hits_x[ihit]-clusters_TM_X[iclus],_hits_y[ihit]-clusters_TM_Y[iclus]);
            // check eta difference
            // TVector3 trackvec(_track_px[itrk],_track_py[itrk],_track_pz[itrk]);
            // trackvec*=(zHC/trackvec.Z());
            // float tracketa = trackvec.Eta();
            // if(!b_TMstudies_x_y_track_filled[icalo])h_TMstudies_x_y_track[icalo]->Fill(trackvec.X(),trackvec.Y());

            // h_TMstudies_dx_vsClsEta_track[0][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[iclus],clusters_TM_Eta[iclus]);
            // h_TMstudies_dx_vsTrkEta_track[0][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[iclus],tracketa);
            // if(fillcurrentlargest){
            //   h_TMstudies_dx_dy_track[1][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[largestClusterID],trackvec.Y()-clusters_TM_Y[largestClusterID]);
            //   h_TMstudies_dx_vsClsEta_track[1][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[largestClusterID],clusters_TM_Eta[largestClusterID]);
            //   h_TMstudies_dx_vsTrkEta_track[1][icalo][ialgo]->Fill(trackvec.X()-clusters_TM_X[largestClusterID],tracketa);
            // }
            if(abs(_hits_x[ihit]-clusters_TM_X[iclus]) < matchingwdw){
              if(abs(_hits_y[ihit]-clusters_TM_Y[iclus]) < matchingwdw){
                ismatched_hit = true;
              }
            }
            // if(abs(trackvec.X()-clusters_TM_X[largestClusterID]) < matchingwdw){
            //   if(abs(trackvec.Y()-clusters_TM_Y[largestClusterID]) < matchingwdw){
            //     ismatched_largest_track = true;
            //   }
            // }
          }
        }

        if(icalo==kFHCAL){
          for(Int_t iclus2=0; iclus2<_nclusters_MA_FEMC; iclus2++){
            if(clusters_TM_X[iclus]==0 && clusters_TM_Y[iclus]==0) continue;
            if(_clusters_MA_FEMC_X[iclus2]==0 && _clusters_MA_FEMC_Y[iclus2]==0) continue;
            // if(clusters_TM_trueID[iclus]==_clusters_MA_FEMC_trueID[iclus2]){
            // if(_mcpart_PDG[clusters_TM_trueID[iclus]]==_mcpart_PDG[_clusters_MA_FEMC_trueID[iclus2]] && _mcpart_PDG[clusters_TM_trueID[iclus]]==211){
              if(abs(_clusters_MA_FEMC_X[iclus2]-clusters_TM_X[iclus]) < 15){
                if(abs(_clusters_MA_FEMC_Y[iclus2]-clusters_TM_Y[iclus]) < 15){

                  h_TMstudies_dx_dy_calomatch[0][icalo][ialgo]->Fill(_clusters_MA_FEMC_X[iclus2]-clusters_TM_X[iclus],_clusters_MA_FEMC_Y[iclus2]-clusters_TM_Y[iclus]);
                  h_TMstudies_E1_E2_calomatch[0][icalo][ialgo]->Fill(_clusters_MA_FEMC_E[iclus2],clusters_TM_E[iclus]);
                  // if(abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==abs(_mcpart_PDG[_clusters_MA_FEMC_trueID[iclus2]]) && abs(_mcpart_PDG[clusters_TM_trueID[iclus]])==211)
                  h_TMstudies_clusSpec_MatchedCalo[icalo]->Fill(clusters_TM_E[iclus]);
                  if(_clusters_MA_FEMC_E[iclus2]>clusters_TM_E[iclus])h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]->Fill(clusters_TM_E[iclus]);
                  if(_clusters_MA_FEMC_E[iclus2]<clusters_TM_E[iclus])h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]->Fill(clusters_TM_E[iclus]);
                }
              }
            // }
          }
        }
        // b_TMstudies_x_y_track_filled[icalo] = true;

        if(ismatched_track){
          h_TMstudies_isMatched[0][icalo]->Fill(3*ialgo+1);
          h_TMstudies_isMatched_E[0][icalo]->Fill(3*ialgo+1,clusters_TM_E[iclus] );
          h_TMstudies_isMatched_PDG[0][icalo]->Fill(3*ialgo+1,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[0][icalo]->Fill(3*ialgo+1);
            h_TMstudies_isMatched_findable_E[0][icalo]->Fill(3*ialgo+1,clusters_TM_E[iclus] );
            h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(3*ialgo+1,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
          }
        }
        if(fillcurrentlargest && ismatched_largest_track){
          h_TMstudies_isMatched[1][icalo]->Fill(3*ialgo+1);
          h_TMstudies_isMatched_E[1][icalo]->Fill(3*ialgo+1,clusters_TM_E[largestClusterID] );
          h_TMstudies_isMatched_PDG[1][icalo]->Fill(3*ialgo+1,abs(_mcpart_PDG[clusters_TM_trueID[largestClusterID]]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[1][icalo]->Fill(3*ialgo+1);
            h_TMstudies_isMatched_findable_E[1][icalo]->Fill(3*ialgo+1,clusters_TM_E[largestClusterID] );
            h_TMstudies_isMatched_findable_PDG[1][icalo]->Fill(3*ialgo+1,abs(_mcpart_PDG[clusters_TM_trueID[largestClusterID]]) );
          }
        }
        if(ismatched_projection){
          h_TMstudies_isMatched[0][icalo]->Fill(3*ialgo+2);
          h_TMstudies_isMatched_E[0][icalo]->Fill(3*ialgo+2,clusters_TM_E[iclus] );
          h_TMstudies_isMatched_PDG[0][icalo]->Fill(3*ialgo+2,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[0][icalo]->Fill(3*ialgo+2);
            h_TMstudies_isMatched_findable_E[0][icalo]->Fill(3*ialgo+2,clusters_TM_E[iclus] );
            h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(3*ialgo+2,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
          }
        }
        if(fillcurrentlargest && ismatched_largest_projection){
          h_TMstudies_isMatched[1][icalo]->Fill(3*ialgo+2);
          h_TMstudies_isMatched_E[1][icalo]->Fill(3*ialgo+2,clusters_TM_E[largestClusterID] );
          h_TMstudies_isMatched_PDG[1][icalo]->Fill(3*ialgo+2,abs(_mcpart_PDG[clusters_TM_trueID[largestClusterID]]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[1][icalo]->Fill(3*ialgo+2);
            h_TMstudies_isMatched_findable_E[1][icalo]->Fill(3*ialgo+2,clusters_TM_E[largestClusterID] );
            h_TMstudies_isMatched_findable_PDG[1][icalo]->Fill(3*ialgo+2,abs(_mcpart_PDG[clusters_TM_trueID[largestClusterID]]) );
          }
        }
        if(ismatched_hit){
          h_TMstudies_isMatched[0][icalo]->Fill(3*ialgo+3);
          h_TMstudies_isMatched_E[0][icalo]->Fill(3*ialgo+3,clusters_TM_E[iclus] );
          h_TMstudies_isMatched_PDG[0][icalo]->Fill(3*ialgo+3,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
          if(isfindable){
            h_TMstudies_isMatched_findable[0][icalo]->Fill(3*ialgo+3);
            h_TMstudies_isMatched_findable_E[0][icalo]->Fill(3*ialgo+3,clusters_TM_E[iclus] );
            h_TMstudies_isMatched_findable_PDG[0][icalo]->Fill(3*ialgo+3,abs(_mcpart_PDG[clusters_TM_trueID[iclus]]) );
          }
        }
      }
    }
  }
  return false;
}

void trackmatchingstudiesSave(){
    // cout << "FHCAL " << nelecinaccFHCAL_TMstudies << " electrons and " << npioninaccFHCAL_TMstudies << " pions in acceptance, and in total " << partinaccFHCAL_TMstudies << " in acc" << endl;
    // cout << "FEMC " << nelecinaccFEMC_TMstudies << " electrons and " << npioninaccFEMC_TMstudies << " pions in acceptance, and in total " << partinaccFEMC_TMstudies << " in acc" << endl;

  // define output file
  TFile* fileOutput = new TFile(Form("%s/output_TMSTUD.root",outputDir.Data()),"RECREATE");
  // fileOutput->mkdir("etabins");
  for(int icalo=0;icalo<_active_calo;icalo++){
    fileOutput->mkdir(Form("%s",str_calorimeter[icalo].Data()));
    fileOutput->cd(Form("%s",str_calorimeter[icalo].Data()));
    for(int icase=0;icase<diffCases_TMStudies;icase++){
      // fileOutput->cd();
      if(h_TMstudies_x_y_track[icalo]&&icase==0)h_TMstudies_x_y_track[icalo]->Write();
      if(h_TMstudies_x_y_projection[icalo]&&icase==0)h_TMstudies_x_y_projection[icalo]->Write();
      if(h_TMstudies_x_y_true[icalo]&&icase==0)h_TMstudies_x_y_true[icalo]->Write();
      if(h_TMstudies_clusSpec[icalo]&&icase==0)h_TMstudies_clusSpec[icalo]->Write();
      if(h_TMstudies_clusSpec_Matched[icalo]&&icase==0)h_TMstudies_clusSpec_Matched[icalo]->Write();
      if(h_TMstudies_clusSpec_MatchedCalo[icalo]&&icase==0)h_TMstudies_clusSpec_MatchedCalo[icalo]->Write();
      if(h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]&&icase==0)h_TMstudies_clusSpec_MatchedCalo_EDiff1[icalo]->Write();
      if(h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]&&icase==0)h_TMstudies_clusSpec_MatchedCalo_EDiff2[icalo]->Write();
      if(h_TMstudies_clusSpec_special[icalo]&&icase==0)h_TMstudies_clusSpec_special[icalo]->Write();
      if(h_TMstudies_clusSpec_special_Matched[icalo]&&icase==0)h_TMstudies_clusSpec_special_Matched[icalo]->Write();
      if(h_TMstudies_isMatched[icase][icalo])h_TMstudies_isMatched[icase][icalo]->Write();
      if(h_TMstudies_isMatched_E[icase][icalo])h_TMstudies_isMatched_E[icase][icalo]->Write();
      if(h_TMstudies_isMatched_PDG[icase][icalo])h_TMstudies_isMatched_PDG[icase][icalo]->Write();
      if(h_TMstudies_isMatched_findable[icase][icalo])h_TMstudies_isMatched_findable[icase][icalo]->Write();
      if(h_TMstudies_isMatched_findable_E[icase][icalo])h_TMstudies_isMatched_findable_E[icase][icalo]->Write();
      if(h_TMstudies_isMatched_findable_PDG[icase][icalo])h_TMstudies_isMatched_findable_PDG[icase][icalo]->Write();
      for(int ialgo=0;ialgo<_active_algo;ialgo++){
        if(h_TMstudies_x_y_clus[icalo][ialgo]&&icase==0) h_TMstudies_x_y_clus[icalo][ialgo]->Write();
        if(h_TMstudies_dx_dy_track[icase][icalo][ialgo]) h_TMstudies_dx_dy_track[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_dy_projection[icase][icalo][ialgo]) h_TMstudies_dx_dy_projection[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_dy_TTLhits[icase][icalo][ialgo]) h_TMstudies_dx_dy_TTLhits[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_dy_calomatch[icase][icalo][ialgo]) h_TMstudies_dx_dy_calomatch[icase][icalo][ialgo]->Write();
        if(h_TMstudies_E1_E2_calomatch[icase][icalo][ialgo]) h_TMstudies_E1_E2_calomatch[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo]) h_TMstudies_dx_vsClsEta_track[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo]) h_TMstudies_dx_vsClsEta_projection[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo]) h_TMstudies_dx_vsTrkEta_track[icase][icalo][ialgo]->Write();
        if(h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo]) h_TMstudies_dx_vsTrkEta_projection[icase][icalo][ialgo]->Write();
      }
      for (Int_t et = 0; et<nEta+1; et++){
        if(h_TMstudies_clusSpec_Eta[icalo][et]&&icase==0)h_TMstudies_clusSpec_Eta[icalo][et]->Write();
        if(h_TMstudies_clusSpec_Matched_Eta[icalo][et]&&icase==0)h_TMstudies_clusSpec_Matched_Eta[icalo][et]->Write();
      }
      // fileOutput->cd("..");
    }
  }
  // write output file
  fileOutput->Write();
  fileOutput->Close();
}
