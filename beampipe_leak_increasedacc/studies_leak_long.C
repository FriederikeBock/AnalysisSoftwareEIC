#include "../common/plottingheader.h"
#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0]))

void fill_histogram(TH1F * const h, TTree * const t, const Float_t true_energy,
		    const Float_t min_value,const Float_t min_towers, const bool normalize, Float_t etalow, Float_t etahigh);
void split_canvas(TCanvas* &cPNG, TString canvName, Int_t numInputs);

void studies_leak_long(
  Double_t jetradiusplot = 0.7,
    TString suffix            = "pdf"
){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem   	              = "e+p: 10#times250 GeV^{2}";
  TString dateForOutput                       = ReturnDateStringForOutput();

  TString outputDir 						                  = Form("plotsleaklong/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  TString inputFolderNames[2] = {"Fun4All_G4_FullDetector_def","Fun4All_G4_FullDetector_HClong"};
  TString namePlotInput[2] = {"default","2m FHCAL"};
  Color_t colorInput[2] = {kGreen+2,kOrange+2};
  Marker_t markerInput[2] = {27,28};


  //************************** Read data **************************************************
  const Int_t nEne = 10;
  // const static Double_t partE[]   = {3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 7.5, 8.0, 9.0,  10., 11., 13., 15., 20., 25., 30., 40., 50., 60., 70., 80., 150.};
  const static Double_t partE[]   = {3.0, 4.0, 5.0, 7.0,  10., 15., 20., 30., 50., 80.};
  Int_t useE[nEne][2]                = {{ 0}};
  TH1F* h_pion_fhcal_E[nEne][2] 	    = {NULL};
  TH2F* h2D_pion_fhcal_dEtadPhi[nEne][2] 	    = {NULL};
  TH1F* h_pion_fhcal_dEta[nEne][2] 	    = {NULL};
  TH1F* h_pion_fhcal_dPhi[nEne][2] 	    = {NULL};
  Int_t activeE[2] = {0};

  Double_t etalowdef = 2.5;
  Double_t etahighdef = 3.0;

  for(Int_t iInp=0; iInp<2;iInp++){
    cout << inputFolderNames[iInp].Data() << endl;
    for(Int_t epart=0; epart<NELEMS(partE); epart++){
      // if(partE[epart]!=15) continue;
      cout << partE[epart] << endl;
      Int_t nevtgen = 2000;
      TChain* t_pion_fhcal_clusters = new TChain("ntp_gshower", "ntp_gshower");
      for(Int_t ii=1; ii<4;ii++){
        if(partE[epart]>70) nevtgen = 1000;
        // if(!gSystem->AccessPathName(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/%d_%d_%1.1f/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),nevtgen,ii,partE[epart]))){
          t_pion_fhcal_clusters->Add(Form("/media/nschmidt/external2/EICsimulationOutputs/%s/%d_%1.1f_%d/G4EICDetector_g4fhcal_eval.root",inputFolderNames[iInp].Data(),nevtgen,partE[epart],ii));
        // }
      }
      if(t_pion_fhcal_clusters){
        useE[epart][iInp] = 1;
        activeE[iInp]++;
        h_pion_fhcal_E[epart][iInp] 	= new TH1F(Form("h_pion_fhcal_E_%d",epart), "", 100, 0.0, 2.0);

        fill_histogram(h_pion_fhcal_E[epart][iInp], t_pion_fhcal_clusters, partE[epart], 0.3,2, true, etalowdef, etahighdef);
        // fill_histogram(h_pion_fhcal_E[epart][iInp], t_pion_fhcal_clusters, partE[epart], 0.3,2, true, etalowdef, etahighdef);
      }
      // h_pion_fhcal_E[epart][iInp]->Scale(1 / h_pion_fhcal_E[epart][iInp]->GetEntries());
    }
  }


  TCanvas* cPNG; // = new TCanvas("cPNG", "", gStyle->GetCanvasDefW() * 3,
	split_canvas(cPNG, "cPNG1", activeE[0]);
  Double_t textSizeSinglePad                    = 0.05;
  Double_t textSizeLabelsPixel                    = 35;
  Double_t textSizeLabelsRel                    = 58./1300;
  TH2F * histEPlotDummy                           = new TH2F("histEPlotDummy","histEPlotDummy",1000,0, 1.5,1000,0, 0.070);
  SetStyleHistoTH2ForGraphs(histEPlotDummy, "#it{E}_{rec} / #it{E}_{true}","norm. counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
  histEPlotDummy->GetXaxis()->SetNoExponent();
  histEPlotDummy->GetYaxis()->SetNdivisions(505,kTRUE);
  histEPlotDummy->GetXaxis()->SetMoreLogLabels(kTRUE);

  TF1 *fit_reso[nEne][2] = {NULL};

  TH1F*  h_pion_resolution_fhcal_E 	= new TH1F("h_pion_resolution_fhcal_E", "", 400, 0.0, 200.0);
  TH1F*  h_pion_resolution_fhcal_1oE 	= new TH1F("h_pion_resolution_fhcal_1oE", "", 500, 0.0, 1.0);
  Int_t padnum =0;
  for(Int_t epart=0; epart<NELEMS(partE); epart++){
    if(!useE[epart][0]) continue;
  	cPNG->SetRightMargin(0.01);
  	cPNG->SetTopMargin(0.01);
  	cPNG->SetLeftMargin(0.08);
  	cPNG->SetBottomMargin(0.1);
  	cPNG->cd(padnum+1);
    histEPlotDummy->DrawCopy();
    drawLatexAdd(Form("%1.1f GeV #pi^{-} in %1.1f<#eta<%1.1f",partE[epart],etalowdef,etahighdef),0.15,0.15,1.3*textSizeLabelsRel,kFALSE,kFALSE,kFALSE);
    drawLatexAdd(Form("rec. w/ FHCAL"),0.95,0.15,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    TLegend* legendHE_gamma           = GetAndSetLegend2(0.30, 0.95-(5*textSizeLabelsRel), 0.95, 0.95,0.7*textSizeLabelsPixel, 1, "", 43, 0.15);
    for(Int_t iInp=0; iInp<2;iInp++){
      if(!useE[epart][iInp]) continue;

      DrawGammaSetMarker(h_pion_fhcal_E[epart][iInp], markerInput[iInp], 2, colorInput[iInp], colorInput[iInp]);
      h_pion_fhcal_E[epart][iInp]->Draw("same");

      fit_reso[epart][iInp] = new TF1(Form("fit_reso_%d",epart), "gaus", 0.4,0.9);//iInp<2 ? fitrange_low[epart]-0.1 : fitrange_low[epart],iInp<2 ? fitrange_hig[epart]-0.1 : fitrange_hig[epart]);

      h_pion_fhcal_E[epart][iInp]->Fit(fit_reso[epart][iInp],"L0RMEQ");
      fit_reso[epart][iInp]->SetNpx(10000);
      DrawGammaSetMarkerTF1(fit_reso[epart][iInp], 1, 2, colorInput[iInp]+2);
      fit_reso[epart][iInp]->Draw("same");
      legendHE_gamma->AddEntry(h_pion_fhcal_E[epart][iInp],Form("%s: mean: %1.2f, #sigma/#it{E}: %1.2f%%",namePlotInput[iInp].Data(),fit_reso[epart][iInp]->GetParameter(1),100*fit_reso[epart][iInp]->GetParameter(2)/fit_reso[epart][iInp]->GetParameter(1)),"pl");
    }
    legendHE_gamma->Draw();

    // drawLatexAdd(Form("mean: %1.2f ",fit_reso[epart]->GetParameter(1)),0.90,0.80,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("sigma: %1.2f ",fit_reso[epart]->GetParameter(2)),0.90,0.75,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd(Form("#sigma/#it{E}: %1.2f %%",100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1)),0.90,0.70,1.3*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    // h_pion_resolution_fhcal_E->SetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(2)/fit_reso[epart]->GetParameter(2),2)+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)));
    // h_pion_resolution_fhcal_E->SetBinError(h_pion_resolution_fhcal_E->FindBin(partE[epart]),h_pion_resolution_fhcal_E->GetBinContent(h_pion_resolution_fhcal_E->FindBin(partE[epart]))*TMath::Sqrt(TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2))); //+TMath::Power(fit_reso[epart]->GetParError(1)/fit_reso[epart]->GetParameter(1),2)
    // h_pion_resolution_fhcal_1oE->SetBinContent(h_pion_resolution_fhcal_1oE->FindBin(1./TMath::Sqrt(partE[epart])),100*fit_reso[epart]->GetParameter(2)/fit_reso[epart]->GetParameter(1));
    padnum++;
 }
  cPNG->Print(Form("%s/HCAL_EtrueFrac_E.%s", outputDir.Data(), suffix.Data()));


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
