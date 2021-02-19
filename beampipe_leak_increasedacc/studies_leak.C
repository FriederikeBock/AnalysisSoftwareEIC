#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value,const Float_t min_towers, const bool normalize, Float_t etalow, Float_t etahigh);
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);

void studies_leak(
  Double_t jetradiusplot = 0.7,
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "e+p: 10#times250 GeV^{2}";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TString inputFolderNames = "Fun4All_G4_FullEtaLeakTestDetector_Pythia"; //"Fun4All_G4_FullDetector_Pythia_nobeampipe","Fun4All_G4_FullDetector_Pythia_beampipe","Fun4All_G4_HighAccDetector_Pythia",};


  //************************** Read data **************************************************
  const Int_t nEne = 22;
  const static Double_t partE[]   = {3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0,  10., 11., 13., 15., 20., 25., 30., 40., 50., 60., 70., 80.};
  Int_t useE[nEne]                = { 0};
  Int_t activeE = 0;


  TH2F* h2D_leakfrac_E_vs_Eta[nEne] 	    = {NULL};
  TH2F* h2D_leakTruefrac_E_vs_Eta[nEne] 	    = {NULL};
  TH1F* h_leakTruefrac_Mean_Eta[nEne] 	    = {NULL};
  TH1F* h_leakfrac_Mean_Eta[nEne] 	    = {NULL};
  // TH2F* h_EtaPhiMap[nEne] 	    = {NULL};
  // TH2F* h2D_NTowers_ECAL_pi[nEne] 	    = {NULL};

// /*
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    cout << partE[epart] << endl;
    // TTree *const t_clusters_fhcal_fulleta =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),partE[epart]), "READ"))->Get("ntp_cluster"); //ntp_cluster
    // TTree *const t_clusters_fhcal_fulleta_iter =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),partE[epart]), "READ"))->Get("ntp_cluster"); //ntp_cluster
    // TTree *const t_towers_fhcal_fulleta =  (TTree *) (new TFile(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),partE[epart]), "READ"))->Get("ntp_tower"); //ntp_cluster

    TChain* t_clusters_fhcal_fulleta = new TChain("ntp_gshower", "ntp_gshower");
    TChain* t_clusters_fhcal_fulleta_iter = new TChain("ntp_gshower", "ntp_gshower");
    TChain* t_towers_fhcal_fulleta = new TChain("ntp_tower", "ntp_tower");
    // t_clusters_fhcal_fulleta->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),partE[epart]));
    // t_clusters_fhcal_fulleta_iter->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),partE[epart]));
    // t_towers_fhcal_fulleta->Add(Form("/media/nschmidt/local/EIC_running/cluster_output/%s/800_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),partE[epart]));
    for(Int_t ii=2; ii<3;ii++){
        t_clusters_fhcal_fulleta->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),ii,partE[epart]));
        t_clusters_fhcal_fulleta_iter->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),ii,partE[epart]));
        t_towers_fhcal_fulleta->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/10000_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames.Data(),ii,partE[epart]));
    }

    // h_EtaPhiMap[epart] 	= new TH2F(Form("h_EtaPhiMap%d",epart), "", 55, -0.5, 54.5, 55, -0.5, 54.5);
    h_leakTruefrac_Mean_Eta[epart] 	= new TH1F(Form("h_leakTruefrac_Mean_Eta%d",epart), "", 60, 2.5, 4.5);
    h_leakfrac_Mean_Eta[epart] 	= new TH1F(Form("h_leakfrac_Mean_Eta%d",epart), "", 60, 2.5, 4.5);
    h2D_leakfrac_E_vs_Eta[epart] 	= new TH2F(Form("h2D_leakfrac_E_vs_Eta_%d",epart), "", 60, 2.5, 4.5,120,0,1.2);
    h2D_leakTruefrac_E_vs_Eta[epart] 	= new TH2F(Form("h2D_leakTruefrac_E_vs_Eta_%d",epart), "", 60, 2.5, 4.5,120,0,1.2);
    // if(partE[epart]<10 || partE[epart]>10) continue;
    if(partE[epart]!=40) continue;
    if(t_clusters_fhcal_fulleta && t_towers_fhcal_fulleta){
        useE[epart] = 1;
        activeE++;

      Float_t clusterEvent,clusterEvent2,towerEvent;
      t_clusters_fhcal_fulleta->SetBranchAddress("event", &clusterEvent);
      t_clusters_fhcal_fulleta_iter->SetBranchAddress("event", &clusterEvent2);
      t_towers_fhcal_fulleta->SetBranchAddress("event", &towerEvent);
      Float_t clusterID,clusterID2;
      t_clusters_fhcal_fulleta->SetBranchAddress("clusterID", &clusterID);
      t_clusters_fhcal_fulleta_iter->SetBranchAddress("clusterID", &clusterID2);
      Float_t clusterE,clusterE2,towerE;
      t_clusters_fhcal_fulleta->SetBranchAddress("e", &clusterE);
      t_clusters_fhcal_fulleta_iter->SetBranchAddress("e", &clusterE2);
      t_towers_fhcal_fulleta->SetBranchAddress("e", &towerE);
      Float_t clusterTrueEta,towerTrueEta;
      t_clusters_fhcal_fulleta->SetBranchAddress("geta", &clusterTrueEta);
      t_towers_fhcal_fulleta->SetBranchAddress("geta", &towerTrueEta);
      Float_t clusterTruePhi,towerTruePhi;
      t_clusters_fhcal_fulleta->SetBranchAddress("gphi", &clusterTruePhi);
      t_towers_fhcal_fulleta->SetBranchAddress("gphi", &towerTruePhi);
      Float_t clusterEta,towerEta;
      t_clusters_fhcal_fulleta->SetBranchAddress("eta", &clusterEta);
      t_towers_fhcal_fulleta->SetBranchAddress("eta", &towerEta);
      Float_t towerIEta,towerIPhi;
      t_towers_fhcal_fulleta->SetBranchAddress("ieta", &towerIEta);
      t_towers_fhcal_fulleta->SetBranchAddress("iphi", &towerIPhi);
      Float_t clusterPhi,towerPhi;
      t_clusters_fhcal_fulleta->SetBranchAddress("phi", &clusterPhi);
      t_towers_fhcal_fulleta->SetBranchAddress("phi", &towerPhi);
      Float_t nTowersCluster;
      t_clusters_fhcal_fulleta->SetBranchAddress("ntowers", &nTowersCluster);

      Float_t lastEvt = -1;
      Double_t maxEcluster = -1;
      Double_t idMaxEcluster = -1;
      Double_t energyLeak = -1;
      Double_t lastEntryUsed = 0;
      for (Int_t i = 0; i < Int_t(t_clusters_fhcal_fulleta->GetEntries()); ++i) {
        // if(i>50) break;
        if (t_clusters_fhcal_fulleta->LoadTree(i) < 0)
          break;

        t_clusters_fhcal_fulleta->GetEntry(i);
        // if(clusterEvent != lastEvt){
        //   lastEvt = clusterEvent;
        //   maxEcluster = -1;
        //   idMaxEcluster = -1;
        //   for (Int_t j = 0; j < Int_t(t_clusters_fhcal_fulleta_iter->GetEntries()); ++j) {
        //     t_clusters_fhcal_fulleta_iter->GetEntry(j);
        //     if(clusterEvent == clusterEvent2){
        //       if(clusterE2>maxEcluster){
        //         maxEcluster = clusterE2;
        //         idMaxEcluster = clusterID2;
        //       }
        //     }
        //   }
        // }
        // if(clusterID !=idMaxEcluster) continue;
        // if(clusterEta<3.0) continue;
        // if(clusterTrueEta<1.0) continue;
        // if(clusterTrueEta>4.5) continue;
        // if(clusterEta>4.0) continue;
        // if(clusterEta<3.5) continue;
        // cout << clusterEvent << endl;
        // cout << "max: " << maxEcluster << "\tID: " << idMaxEcluster << endl;
        // cout << "cur: " << clusterE << "\tID: " << clusterID << "\teta: " << clusterEta<< endl;
        // cout << clusterEta << endl;
        // Loop over high eta cells to get leak energy
        energyLeak = 0;
        for (Int_t k = lastEntryUsed; k < Int_t(t_towers_fhcal_fulleta->GetEntries()); ++k) {
          t_towers_fhcal_fulleta->GetEntry(k);
          // h_EtaPhiMap[epart]->Fill(towerIEta,towerIPhi);
          if(towerEvent == clusterEvent){
            if(towerIEta>25 && towerIEta<29 && towerIPhi>25 && towerIPhi<29){
              energyLeak += towerE;
              lastEntryUsed = k;
            }
          } else if(towerEvent>clusterEvent){
            break;
          }
        }
        h2D_leakfrac_E_vs_Eta[epart]->Fill(clusterEta,energyLeak/clusterE);
        h2D_leakTruefrac_E_vs_Eta[epart]->Fill(clusterTrueEta,energyLeak/partE[epart]);
        // cout << "leak: " << energyLeak << "\tfrac: " << energyLeak/clusterE << endl;
        // (measured_energy1 > 0.1) // above MIP
        // && (measured_energy1 > 0.1) // above MIP
        // && nTowers1>=2 && nTowers2>=2 // actual cluster with >1 tower
        // ){

        // }
      }
    }
    for (Int_t i=1; i < h2D_leakTruefrac_E_vs_Eta[epart]->GetNbinsX(); i++){
      TH1D* projectionYdummy = (TH1D*)h2D_leakTruefrac_E_vs_Eta[epart]->ProjectionY(Form("projectionYdummy%d%d",i,epart), i,i+1,"e");
      h_leakTruefrac_Mean_Eta[epart]->SetBinContent(i,projectionYdummy->GetMean());
      h_leakTruefrac_Mean_Eta[epart]->SetBinError(i,projectionYdummy->GetMeanError());
    }
    for (Int_t i=1; i < h2D_leakfrac_E_vs_Eta[epart]->GetNbinsX(); i++){
      TH1D* projectionYdummy = (TH1D*)h2D_leakfrac_E_vs_Eta[epart]->ProjectionY(Form("projectionYdummy%d%d",i,epart), i,i+1,"e");
      h_leakfrac_Mean_Eta[epart]->SetBinContent(i,projectionYdummy->GetMean());
      h_leakfrac_Mean_Eta[epart]->SetBinError(i,projectionYdummy->GetMeanError());
    }
  }

  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
  Int_t padnum =0;

// 	split_canvas(cPNG, "cPNG1", activeE);
//   TH2F* histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy", 55, -0.5, 54.5, 55, -0.5, 54.5);
//   SetStyleHistoTH2ForGraphs(histEPlotDummy, "#eta_{tower}","#phi_{tower}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
//   histEPlotDummy->GetXaxis()->SetNoExponent();
//   histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
//   histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

//   padnum =0;
//   for(Int_t epart=0; epart<NELEMS(partE); epart++){
//     if(!useE[epart]) continue;
//   	cPNG->cd(padnum+1);
//     histEPlotDummy->DrawCopy();
//     h_EtaPhiMap[epart]->Draw("col");
//     drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.90,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

//     padnum++;
//  }
//   cPNG->Print(Form("%s/EtaPhiMap.%s", outputDir.Data(), suffix.Data()));
  

	split_canvas(cPNG, "cPNG1", activeE);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
  	cPNG->SetLogz();
  	cPNG->SetLeftMargin(0.1);
  	cPNG->SetRightMargin(0.1);
  	cPNG->SetBottomMargin(0.1);
  	cPNG->cd(padnum+1);
    // DrawGammaSetMarker(h2D_leakfrac_E_vs_Eta[epart], markerInput[0], 2, colorInput[0], colorInput[0]);
    SetStyleHistoTH2ForGraphs(h2D_leakfrac_E_vs_Eta[epart], "#eta_{cluster}","#it{E}_{leak} / #it{E}_{cluster}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    h2D_leakfrac_E_vs_Eta[epart]->GetYaxis()->SetRangeUser(0,1.2);
    h2D_leakfrac_E_vs_Eta[epart]->Draw("colz");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.85,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("#it{E}_{leak} = #sum_{#eta>4} #it{E}_{tower}",0.15,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd("#it{E}_{leak} = #it{E}^{#eta>4}_{shower}",0.15,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    DrawGammaSetMarker(h_leakfrac_Mean_Eta[epart], 20, 1, kOrange+2, kOrange+2);
    h_leakfrac_Mean_Eta[epart]->Draw("same");
    padnum++;
 }
  cPNG->Print(Form("%s/ELeakFracEta_E.%s", outputDir.Data(), suffix.Data()));

	split_canvas(cPNG, "cPNG1", activeE);
  padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart]) continue;
    TCanvas* cECCE = new TCanvas("cECCE","",0,0,1100,1000);
    DrawGammaCanvasSettings( cECCE, 0.11, 0.01, 0.01, 0.105);
  	// cPNG->SetLeftMargin(0.1);
  	// cPNG->SetRightMargin(0.02);
  	// cPNG->SetBottomMargin(0.11);
  	// cPNG->cd(padnum+1);
    SetStyleHistoTH2ForGraphs(h2D_leakTruefrac_E_vs_Eta[epart], "#eta_{#pi}^{true}","#it{E}^{rec}_{leak} / #it{E}^{#pi}_{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.99);
    // DrawGammaSetMarker(h2D_leakfrac_E_vs_Eta[epart], markerInput[0], 2, colorInput[0], colorInput[0]);
    h2D_leakTruefrac_E_vs_Eta[epart]->GetYaxis()->SetRangeUser(0,0.79);
    h2D_leakTruefrac_E_vs_Eta[epart]->GetXaxis()->SetRangeUser(3.01,4.0);
    h2D_leakTruefrac_E_vs_Eta[epart]->Draw("col");
    drawLatexAdd(Form("%1.1f GeV #pi^{-} on FHCAL",partE[epart]),0.15,0.90,1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    // drawLatexAdd("#it{E}_{leak} = #sum_{#eta>4} #it{E}_{tower}",0.15,0.85,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd("#it{E}_{leak} = #it{E}^{#eta>4}_{shower}",0.15,0.83,1*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    DrawGammaSetMarker(h_leakTruefrac_Mean_Eta[epart], 20, 3, kOrange+2, kOrange+2);
    h_leakTruefrac_Mean_Eta[epart]->Draw("same");
    TLegend* legendJES  = GetAndSetLegend2(0.13, 0.83-(2*textSizeLabelsRel), 0.6, 0.83-(1*textSizeLabelsRel),textSizeLabelsPixel, 1, "", 43, 0.15);
    legendJES->AddEntry(h_leakTruefrac_Mean_Eta[epart],"mean transverse leakage","p");
    legendJES->Draw();
    padnum++;
    cECCE->Print(Form("%s/ELeakTrueFracTrueEta_E.%s", outputDir.Data(), suffix.Data()));
 }
  // cPNG->Print(Form("%s/ELeakTrueFracTrueEta_E.%s", outputDir.Data(), suffix.Data()));
  // */
}
          // h2D_NTowers_HCAL_pi[epart][iInp]->Fill(measured_energy1, nTowers1);


void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value, const Float_t min_towers, const bool normalize, Float_t etalow, Float_t etahigh)
{
	Float_t measured_energy;
//      Float_t true_energy;
	t->SetBranchAddress("e", &measured_energy);
//      t->SetBranchAddress("ge", &true_energy);
	Float_t true_eta;
	t->SetBranchAddress("geta", &true_eta);
	Float_t clusterID;
	t->SetBranchAddress("clusterID", &clusterID);
	Float_t nTowers;
	t->SetBranchAddress("ntowers", &nTowers);
	Float_t eventID;
	t->SetBranchAddress("event", &eventID);

	Int_t nentries = Int_t(t->GetEntries());
  Float_t lastEvt = 0;
	for (Int_t i = 0; i < nentries; ++i) {
		if (t->LoadTree(i) < 0)
			break;

		t->GetEntry(i);
    // if(eventID != lastEvt)cout << eventID << "\tE/Etrue: " << measured_energy/true_energy << endl;
    // else cout << "\tE/Etrue: " << measured_energy/true_energy << endl;
		if (
      // ((true_eta > -0.5 && true_eta < 0.5)
      // || (true_eta > -3 && true_eta < -2)
      // || 
      (true_eta > etalow && true_eta < etahigh) //)  //1.7 < 3
    && (measured_energy > min_value && true_energy > 0.1) // above MIP
    // && clusterID==1
    && nTowers>=min_towers // actual cluster with >1 tower
    ){
			h->Fill(measured_energy / true_energy);
    }
    // if(eventID != lastEvt)lastEvt = eventID;

	}
	if (normalize)
		h->Scale(1 / h->GetEntries());

	h->SetXTitle("E_{cluster} / E_{true}");
	h->SetYTitle("entries / #scale[0.5]{#sum} entries      ");
}

void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs)
{
  if(numInputs<2){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 1, gStyle->GetCanvasDefH()*1);
	  // cPNG->Divide(2, 2);
  }else if(numInputs<5){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 2, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(2, 2);
  }else if(numInputs<7){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*2);
	  cPNG->Divide(3, 2);
  }else if(numInputs<10){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 3, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(3, 3);
  } else if(numInputs<13){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*3);
	  cPNG->Divide(4, 3);
  } else if(numInputs<17){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 4, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(4, 4);
  } else if(numInputs<21){
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*4);
	  cPNG->Divide(5, 4);
  } else {
    cPNG = new TCanvas(canvName.Data(), "", gStyle->GetCanvasDefW() * 5, gStyle->GetCanvasDefH()*5);
	  cPNG->Divide(5, 5);
  }

}
